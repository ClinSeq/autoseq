import logging
import os
import uuid

import sys

from bcbio.variation import vcfutils
from pypedream.job import Job, repeat, required, optional, conditional


class Mutect2(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_tumor = None
        self.input_normal = None
        self.reference_sequence = None
        self.target_regions = None
        self.cosmic = None
        self.dbsnp = None
        self.output = None
        self.scratch = "/tmp/"
        self.jobname = "mutect2"

    def command(self):
        if not self.output.endswith("gz"):
            raise ValueError("Output needs to be gzipped: {}".format(self.output))
        tmpf = "{scratch}/mutect2-{uuid}.vcf.gz".format(scratch=self.scratch, uuid=uuid.uuid4())
        mutect_cmd = "mutect2 -R " + self.reference_sequence + \
                     required("--input_file:tumor ", self.input_tumor) + \
                     optional("--input_file:normal ", self.input_normal) + \
                     required("--out ", tmpf) + \
                     optional("-L ", self.target_regions) + \
                     optional("-nct ", self.threads) + \
                     optional("--dbsnp ", self.dbsnp) + \
                     optional("--cosmic ", self.cosmic)
        leftaln_cmd = "gzip -cd {} ".format(tmpf) + \
                      " | " + vt_split_and_leftaln(self.reference_sequence) + \
                      " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        # TODO also filter so only variants with PASS remain
        rmtmp_cmd = "rm {}".format(tmpf)
        return " && ".join([mutect_cmd, leftaln_cmd, rmtmp_cmd])


class Freebayes(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bams = None
        self.tumorid = None
        self.normalid = None
        self.reference_sequence = None
        self.target_bed = None
        self.somatic_only = False
        self.params = "--pooled-discrete --pooled-continuous --genotype-qualities --report-genotype-likelihood-max --allele-balance-priors-off"
        self.min_coverage = 20
        self.min_alt_frac = 0.01
        self.use_harmonic_indel_quals = False
        self.scratch = "/tmp/"
        self.output = ""
        self.jobname = "freebayes-somatic"

    def command(self):
        if not self.tumorid and not self.normalid and self.somatic_only:
            logging.error("Both tumorid and normal id are required for somatic calling.")
            raise ValueError

        regions_file = "{scratch}/{uuid}.regions".format(scratch=self.scratch, uuid=uuid.uuid4())
        bed_to_regions_cmd = "cat {} | bed_to_regions.py > {}".format(self.target_bed, regions_file)
        freebayes_cmd = "freebayes-parallel {} {} ".format(regions_file, self.threads) + \
                        required("-f ", self.reference_sequence) + " --use-mapping-quality " + \
                        optional("--min-alternate-fraction ", self.min_alt_frac) + \
                        optional("--min-coverage ", self.min_coverage) + \
                        conditional(self.use_harmonic_indel_quals, "--harmonic-indel-quality") + \
                        optional("", self.params) + \
                        repeat(" ", self.input_bams) + \
                        """| bcftools filter -i 'ALT="<*>" || QUAL > 5' """ + \
                        "| filter_erroneus_alt.py -V /dev/stdin " + \
                        conditional(self.somatic_only,
                                    " | freebayes-somatic-filter.py -V - " +
                                    " -tumorid {}".format(self.tumorid) +
                                    " -normalid {}".format(self.normalid)) + \
                        " | " + vt_split_and_leftaln(self.reference_sequence) + \
                        " | vcfuniq " + \
                        " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        # reason for 'vcfuniq': freebayes sometimes report duplicate variants that need to be uniqified.
        rm_regions_cmd = "rm {}".format(regions_file)
        return " && ".join([bed_to_regions_cmd, freebayes_cmd, rm_regions_cmd])


class VarDict(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 target_bed=None, output=None, min_alt_frac=0.1):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output = output
        self.min_alt_frac = min_alt_frac

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)

        freq_filter = (" bcftools filter -m '+' -s 'REJECT' -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                       "| %s -c 'import bcbio.variation.vardict; import sys; print bcbio.variation.vardict.depth_freq_filter(sys.stdin.read(), %s, \"%s\")' " %
                       (sys.executable, 0, 'bwa'))

        somatic_filter = (" sed 's/\\.*Somatic\\\"/Somatic/' "  # changes \".*Somatic\" to Somatic
                          "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                          "| %s -c 'import bcbio.variation.freebayes; import sys; print bcbio.variation.freebayes.call_somatic(sys.stdin.read())' " % sys.executable)

        fix_ambig = vcfutils.fix_ambiguous_cl()
        remove_dup = vcfutils.remove_dup_cl()

        cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
              optional("-f ", self.min_alt_frac) + \
              required("-N ", self.tumorid) + \
              " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
              " -c 1 -S 2 -E 3 -g 4 -Q 10 " + required("", self.target_bed) + \
              " | testsomatic.R " + \
              " | var2vcf_paired.pl -P 0.9 -m 4.25 -M " + required("-f ", self.min_alt_frac) + \
              " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
              " | " + freq_filter + " | " + somatic_filter + " | " + fix_ambig + " | " + remove_dup + \
              " | vcfstreamsort -w 1000 " + \
              " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        return cmd


class VEP(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_vcf = None
        self.reference_sequence = None
        self.vep_bin_dir = None
        self.vep_cache_dir = None
        self.jobname = "vep"

    def command(self):
        bgzip = ""
        fork = ""
        if self.threads > 1:  # vep does not accept "--fork 1", so need to check.
            fork = "--fork {}".format(self.threads)
        if self.output_vcf.endswith('gz'):
            bgzip = "| bgzip"
        required("", self.vep_bin_dir)  # assert that vep_bin_dir is set
        return "perl {}/variant_effect_predictor.pl ".format(self.vep_bin_dir) + \
               required('-i ', self.input_vcf) + \
               " --vcf --offline --output_file STDOUT " + \
               required('--dir ', self.vep_cache_dir) + \
               required('--fasta ', self.reference_sequence) + \
               " --check_alleles --check_existing " + \
               " --total_length --allele_number --no_escape --no_stats --everything " + \
               fork + \
               bgzip + " > " + required("", self.output_vcf) + \
               " && tabix -p vcf {}".format(self.output_vcf)


class VcfAddSample(Job):
    """
    Add DP, RO and AO tags for a new sample to a VCF, filter low-qual variants on the fly
    """

    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.input_bam = None
        self.samplename = None
        self.filter_hom = True
        self.output = None
        self.scratch = "/tmp"
        self.jobname = "vcf-add-sample"

    def command(self):
        filt_vcf = "{scratch}/{uuid}.vcf.gz".format(scratch=self.scratch, uuid=uuid.uuid4())
        bgzip = ""
        tabix = ""
        if self.output.endswith('gz'):
            bgzip = "| bgzip"
            tabix = " && tabix -p vcf {}".format(self.output)

        filt_vcf_cmd = "vcf_filter.py --no-filtered " + required("", self.input_vcf) + " sq --site-quality 5 " + \
                       "|bgzip" + " > " + filt_vcf
        vcf_add_sample_cmd = "vcf_add_sample.py " + \
                             conditional(self.filter_hom, "--filter_hom") + \
                             required("--samplename ", self.samplename) + \
                             filt_vcf + " " + \
                             required("", self.input_bam) + \
                             bgzip + " > " + self.output + tabix
        rm_filt_cmd = "rm " + filt_vcf
        return " && ".join([filt_vcf_cmd, vcf_add_sample_cmd, rm_filt_cmd])


def vt_split_and_leftaln(reference_sequence):
    return "vt decompose -s - |vt normalize -r {} -".format(reference_sequence)


class VcfFilter(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.filter = None
        self.output = None
        self.jobname = "vcffilter"

    def command(self):
        return "vcffilter " + \
               required("-f ", self.filter) + \
               required(" ", self.input) + \
               "| bgzip " + required(" > ", self.output) + \
               " && tabix -p vcf {output}".format(output=self.output)


class CurlSplitAndLeftAlign(Job):
    def __init__(self):
        Job.__init__(self)
        self.remote = None
        self.input_reference_sequence = None
        self.input_reference_sequence_fai = None
        self.output = None
        self.jobname = "curl-split-leftaln"

    def command(self):
        required("", self.input_reference_sequence_fai)
        return "curl " + \
               required(" ", self.remote) + \
               "| gzip -d |" + vt_split_and_leftaln(self.input_reference_sequence) + \
               "| bgzip " + required(" > ", self.output) + \
               " && tabix -p vcf {output}".format(output=self.output)


class InstallVep(Job):
    def __init__(self):
        Job.__init__(self)
        self.ensembl_version = None
        self.output_bin_dir = None
        self.output_cache_dir = None
        self.jobname = "install-vep"

    def command(self):
        return "mkdir -p " + self.output_bin_dir + " && " + \
               "mkdir -p " + self.output_cache_dir + " && " + \
               "cd " + self.output_bin_dir + " && " + \
               "curl -LO https://github.com/Ensembl/ensembl-tools/archive/release/" + str(
            self.ensembl_version) + ".tar.gz" + " && " + \
               """ tar -zxf {}.tar.gz  --starting-file variant_effect_predictor --transform='s|.*/|./|g' && """.format(
                   self.ensembl_version) + \
               "export PERL5LIB=$PERL5LIB:" + self.output_bin_dir + " && " + \
               "perl INSTALL.pl --AUTO acf --SPECIES homo_sapiens_vep " + \
               required("--VERSION ", self.ensembl_version) + \
               " --CONVERT " + \
               required("--CACHEDIR ", self.output_cache_dir) + \
               required("--DESTDIR ", self.output_bin_dir)
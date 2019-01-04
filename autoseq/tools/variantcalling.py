import logging
import sys
import uuid

from pypedream.job import Job, repeat, required, optional, conditional
from autoseq.util.clinseq_barcode import *
from autoseq.util.vcfutils import vt_split_and_leftaln, fix_ambiguous_cl, remove_dup_cl

class HaplotypeCaller(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bam = None
        self.reference_sequence = None
        self.java_options = "--java-options -Xmx4g"
        self.output = None
        self.jobname = "gatk-haplotype-somatic"

    def command(self):
        haplotypecaller_cmd = "gatk {} HaplotypeCaller ".format(self.java_options) + \
                        required(" -R ", self.reference_sequence) + \
                        required(" -I ", self.input_bam) + \
                        required(" -O ", self.output) 

        return haplotypecaller_cmd


class VarDict(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 reference_dict=None, target_bed=None, output=None, min_alt_frac=0.1, min_num_reads=None,
                 blacklist_bed=None):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.reference_dict = reference_dict
        self.target_bed = target_bed
        self.blacklist_bed = blacklist_bed
        self.output = output
        self.min_alt_frac = min_alt_frac
        self.min_num_reads = min_num_reads

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)

        freq_filter = (" bcftools filter -e 'STATUS !~ \".*Somatic\"' 2> /dev/null "
                       "| %s -c 'from autoseq.util.bcbio import depth_freq_filter_input_stream; import sys; print depth_freq_filter_input_stream(sys.stdin, %s, \"%s\")' " %
                       (sys.executable, 0, 'bwa'))

        somatic_filter = (" sed 's/\\.*Somatic\\\"/Somatic/' "  # changes \".*Somatic\" to Somatic
                          "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' "
                          "| %s -c 'from autoseq.util.bcbio import call_somatic; import sys; print call_somatic(sys.stdin.read())' " % sys.executable)

        blacklist_filter = " | intersectBed -a . -b {} | ".format(self.blacklist_bed)

        cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
              optional("-f ", self.min_alt_frac) + \
              required("-N ", self.tumorid) + \
              optional("-r ", self.min_num_reads) + \
              " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
              " -c 1 -S 2 -E 3 -g 4 -Q 10 " + required("", self.target_bed) + \
              " | testsomatic.R " + \
              " | var2vcf_paired.pl -P 0.9 -m 4.25 -M " + required("-f ", self.min_alt_frac) + \
              " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
              " | " + freq_filter + " | " + somatic_filter + " | " + fix_ambiguous_cl() + " | " + remove_dup_cl() + \
              " | vcfstreamsort -w 1000 " + \
              " | " + vt_split_and_leftaln(self.reference_sequence) + \
              " | bcftools view --apply-filters .,PASS " + \
              " | vcfsorter.pl {} /dev/stdin ".format(self.reference_dict) + \
              conditional(blacklist_filter, self.blacklist_bed) + \
              " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        return cmd

class StrelkaSomatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 target_bed=None, output_dir=None ):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output_dir = output_dir
        
    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration
        configure_strelkasomatic = "configureStrelkaSomaticWorkflow.py " + \
                                    " --normalBam " + self.input_normal + \
                                    " --tumorBam " + self.input_tumor + \
                                    " --ref " +  self.reference_sequence + \
                                    " --runDir " + self.output_dir
        cmd = configure_strelkasomatic + " && " + self.output_dir+"/runWorkflow.py -m local -j 20"
        return cmd

class StrelkaGermline(Job):
    def __init__(self, input_bam=None, normalid=None, reference_sequence=None,
                 target_bed=None, output_dir=None ):
        Job.__init__(self)
        self.input_bam = input_bam
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output_dir = output_dir
        
    def command(self):
        required("", self.input_bam)
        required("", self.reference_sequence)

        # configuration
        configure_strelkagermline = "configureStrelkaGermlineWorkflow.py " + \
                                    " --bam " + self.input_bam + \
                                    " --ref " +  self.reference_sequence + \
                                    " --runDir " + self.output_dir
        cmd = configure_strelkagermline + " && " + self.output_dir+"/runWorkflow.py -m local -j 20"
        return cmd

class Mutect2Somatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 target_bed=None, output=None, bamout=None, exac=None, tumor_getpileupsummaries_table=None, tumor_calculatecontamination_table=None, output_filtered=None ):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output = output
        self.bamout = bamout
        self.exac_genome_vcf = exac
        self.tumor_getpileupsummaries_table = tumor_getpileupsummaries_table
        self.tumor_calculatecontamination_table = tumor_calculatecontamination_table
        self.output_filtered = output_filtered

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration
        # "-L " + \ We can update Interval List Once confirmed with Rebecka
        # Call somatic short variants and generate a bamout with Mutect2
        mutectsomatic_cmd = " gatk --java-options '-Xmx2g' Mutect2 " + \
                                    " -R " +  self.reference_sequence + \
                                    " -I " + self.input_tumor + \
                                    " -I " + self.input_normal + \
                                    " -tumor " + self.tumorid + \
                                    " -normal " + self.normalid + \
                                    " --af-of-alleles-not-in-resource 0.0000025 " + \
                                    " --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter " + \
                                    " -bamout " + self.bamout + \
                                    " -O " + self.output
        
        # Estimate cross-sample contamination using GetPileupSummaries and CalculateContamination.
        # Run GetPileupSummaries on the tumor BAM to summarize read support for a set number of known variant sites.
        mutect_getpileup_sum = "gatk GetPileupSummaries " + \
                                  "-I " + self.input_tumor + \
                                  "-V " + self.exac_genome_vcf + \
                                  "-O " + self.tumor_getpileupsummaries_table

        # Estimate contamination with CalculateContamination.
        mutect_cal_contamination = "gatk CalculateContamination " + \
                                      "-I " + self.tumor_getpileupsummaries_table + \
                                      "-O " + self.tumor_calculatecontamination_table 

        # Filter for confident somatic calls using FilterMutectCalls 
        filter_mutect_calls = "gatk FilterMutectCalls " + \
                                "-V " + self.output + \
                                "--contamination-table " + self.tumor_calculatecontamination_table + \
                                "-O "  + self.output_filtered

        return " && ".join([mutectsomatic_cmd, mutect_getpileup_sum, mutect_cal_contamination, filter_mutect_calls])


class Varscan2Somatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 target_bed=None, output=None, bamout=None ):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output = output
        
    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration
        # "-L " + \ We can update Interval List Once confirmed with Rebecka
        normal_mpileup_cmd = "samtools mpileup -C50 -f " + self.reference_sequence + " " + self.input_normal + " > " + self.output + "/" + self.normalid +".pileup "
        tumor_mpileup_cmd = "samtools mpileup -C50 -f " + self.reference_sequence + " " + self.input_tumor + " > " + self.output + "/" +  self.tumorid +".pileup "

        varscan_cmd = "java -jar VarScan.v2.4.0.jar somatic " + self.output + "/" + self.normalid +".pileup " + self.output + "/" +  self.tumorid +".pileup " +  \
                      self.output + "/"+ normalid + "-" + tumorid + "-varscan-somatic --output-vcf" 

        return normal_mpileup_cmd + " && " + tumor_mpileup_cmd + " && " + varscan_cmd


class VarDictForPureCN(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 reference_dict=None, target_bed=None, output=None, min_alt_frac=0.1, min_num_reads=None, dbsnp=None):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.reference_dict = reference_dict
        self.target_bed = target_bed
        self.output = output
        self.min_alt_frac = min_alt_frac
        self.min_num_reads = min_num_reads
        self.dbsnp = dbsnp

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)

        tmp_vcf = "{scratch}/{uuid}.vcf.gz".format(scratch=self.scratch, uuid=uuid.uuid4())

        # run vardict without removing non-somatic variants, and adding "SOMATIC" INFO field for somatic variants
        vardict_cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
                      optional("-f ", self.min_alt_frac) + \
                      required("-N ", self.tumorid) + \
                      optional("-r ", self.min_num_reads) + \
                      " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
                      " -c 1 -S 2 -E 3 -g 4 -Q 10 " + required("", self.target_bed) + \
                      " | testsomatic.R " + \
                      " | var2vcf_paired.pl -P 0.9 -m 4.25 " + required("-f ", self.min_alt_frac) + \
                      " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
                      " | " + fix_ambiguous_cl() + " | " + remove_dup_cl() + \
                      " | sed 's/Somatic;/Somatic;SOMATIC;/g' " + \
                      " | sed '/^#CHROM/i ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">' " + \
                      " | vcfstreamsort -w 1000 " + \
                      " | bcftools view --apply-filters .,PASS " + \
                      " | vcfsorter.pl {} /dev/stdin ".format(self.reference_dict) + \
                      " | bgzip > " + tmp_vcf + " && tabix -p vcf " + tmp_vcf

        # annotate variants with dbSNP id
        annotate_cmd = "bcftools annotate --annotation {} --columns ID ".format(self.dbsnp) + \
                       " --output-type z --output {} ".format(self.output) + tmp_vcf + \
                       " && tabix -p vcf {}".format(self.output)

        # remove temporary vcf and tabix
        rm_tmp_cmd = "rm " + tmp_vcf + "*"

        return " && ".join([vardict_cmd, annotate_cmd, rm_tmp_cmd])


class VEP(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_vcf = None
        self.reference_sequence = None
        self.vep_dir = None
        self.jobname = "vep"
        self.additional_options = ""

    def command(self):
        bgzip = ""
        fork = ""
        if self.threads > 1:  # vep does not accept "--fork 1", so need to check.
            fork = " --fork {} ".format(self.threads)
        if self.output_vcf.endswith('gz'):
            bgzip = " | bgzip "

        cmdstr = "variant_effect_predictor.pl --vcf --output_file STDOUT " + \
                 self.additional_options + required("--dir ", self.vep_dir) + \
                 required("--fasta ", self.reference_sequence) + \
                 required("-i ", self.input_vcf) + \
                 " --check_existing  --total_length --allele_number " + \
                 " --no_escape --no_stats --everything --offline " + \
                 fork + bgzip + " > " + required("", self.output_vcf) + \
                 " && tabix -p vcf {}".format(self.output_vcf)

        return cmdstr


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


class VcfFilter(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.filter = None
        self.output = None
        self.jobname = "vcffilter"

    def command(self):
        return "zcat" + \
               required(" ", self.input) + \
               "| vcffilter " + \
               required("-f ", self.filter) + \
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
        return "curl -L " + \
               required(" ", self.remote) + \
               "| gzip -d |" + vt_split_and_leftaln(self.input_reference_sequence, allow_ref_mismatches=True) + \
               "| bgzip " + required(" > ", self.output) + \
               " && tabix -p vcf {output}".format(output=self.output)


class InstallVep(Job):
    def __init__(self):
        Job.__init__(self)
        self.output_dir = None
        self.jobname = "fetch-vep-cache"

    def command(self):
        return "vep_install.pl --SPECIES homo_sapiens_vep --AUTO c --ASSEMBLY GRCh37 --NO_HTSLIB " + \
               required("--CACHEDIR ", self.output_dir) + \
               " && vep_convert_cache.pl " + required("--dir ", self.output_dir) + \
               " --species homo_sapiens --version 83_GRCh37"


def call_somatic_variants(pipeline, cancer_bam, normal_bam, cancer_capture, normal_capture,
                          target_name, outdir, callers=['vardict', 'strelka'],
                          min_alt_frac=0.1, min_num_reads=None):
    """
    Configuring calling of somatic variants on a given pairing of cancer and normal bam files,
    using a set of specified algorithms.

    :param pipeline: The analysis pipeline for which to configure somatic calling.
    :param cancer_bam: Location of the cancer sample bam file
    :param normal_bam: Location of the normal sample bam file
    :param cancer_capture: A UniqueCapture item identifying the cancer sample library capture 
    :param normal_capture: A UniqueCapture item identifying the normal sample library capture
    :param target_name: The name of the capture panel used
    :param outdir: Output location
    :param callers: List of calling algorithms to use - can include 'vardict' and/or 'freebayes'
    :param min_alt_frac: The minimum allelic fraction value in order to retain a called variant 
    :return: A dictionary with somatic caller name as key and corresponding output file location as value
    """
    cancer_capture_str = compose_lib_capture_str(cancer_capture)
    normal_capture_str = compose_lib_capture_str(normal_capture)
    normal_sample_str = compose_sample_str(normal_capture)
    tumor_sample_str = compose_sample_str(cancer_capture)

    d = {}
    if 'freebayes' in callers:
        freebayes = Freebayes()
        freebayes.input_bams = [cancer_bam, normal_bam]
        freebayes.tumorid = cancer_capture_str
        freebayes.normalid = normal_capture_str
        freebayes.somatic_only = True
        freebayes.reference_sequence = pipeline.refdata['reference_genome']
        freebayes.target_bed = pipeline.refdata['targets'][target_name]['targets-bed-slopped20']
        freebayes.threads = pipeline.maxcores
        freebayes.min_alt_frac = min_alt_frac
        freebayes.scratch = pipeline.scratch
        freebayes.jobname = "freebayes-somatic/{}".format(cancer_capture_str)
        freebayes.output = "{}/variants/{}-{}.freebayes-somatic.vcf.gz".format(outdir, cancer_capture_str, normal_capture_str)
        pipeline.add(freebayes)
        d['freebayes'] = freebayes.output

    capture_name = pipeline.get_capture_name(cancer_capture.capture_kit_id)
    blacklist_bed = pipeline.refdata["targets"][capture_name]["blacklist-bed"]

    if 'vardict' in callers:
        vardict = VarDict(input_tumor=cancer_bam, input_normal=normal_bam, tumorid=tumor_sample_str,
                          normalid=normal_sample_str,
                          reference_sequence=pipeline.refdata['reference_genome'],
                          reference_dict=pipeline.refdata['reference_dict'],
                          target_bed=pipeline.refdata['targets'][target_name]['targets-bed-slopped20'],
                          output="{}/variants/{}-{}.vardict-somatic.vcf.gz".format(outdir, cancer_capture_str, normal_capture_str),
                          min_alt_frac=min_alt_frac, min_num_reads=min_num_reads,
                          blacklist_bed=blacklist_bed
                          )

        vardict.jobname = "vardict/{}".format(cancer_capture_str)
        pipeline.add(vardict)
        d['vardict'] = vardict.output


    if 'strelka' in callers:
        strelka_somatic = StrelkaSomatic(input_tumor=cancer_bam, input_normal=normal_bam, tumorid=tumor_sample_str,
                          normalid=normal_sample_str,
                          reference_sequence=pipeline.refdata['reference_genome'],
                          output_dir="{}/variants/{}-{}-strelka-somatic".format(outdir, cancer_capture_str, normal_capture_str)
                          )
        strelka_somatic.jobname = "strelka-somatic-workflow/{}".format(cancer_capture_str)
        pipeline.add(strelka_somatic)
        d['strelka'] = strelka_somatic.output_dir+"/results/variants/somatic.snvs.vcf.gz"

    if 'mutect2' in callers:
        mutect_somatic = Mutect2Somatic(input_tumor=cancer_bam, input_normal=normal_bam, tumorid=tumor_sample_str,
                          normalid=normal_sample_str,
                          reference_sequence=pipeline.refdata['reference_genome'],
                          output="{}/variants/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, cancer_capture_str, normal_capture_str),
                          bamout="{}/bams/{}-{}-mutect.bam".format(outdir, cancer_capture_str, normal_capture_str),
                          exac=pipeline.refdata['exac'],
                          tumor_getpileupsummaries_table= "{}/variants/{}-mutect-tumor-pileupsummary-table".format(outdir, cancer_capture_str),
                          tumor_calculatecontamination_table= "{}/variants/{}-mutect-tumor-contamination-table".format(outdir, cancer_capture_str),
                          output_filtered="{}/variants/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, cancer_capture_str, normal_capture_str)
                          )
        mutect_somatic.jobname = "mutect2-somatic/{}".format(cancer_capture_str)
        pipeline.add(mutect_somatic)
        d['mutect'] = mutect_somatic.output

    return d

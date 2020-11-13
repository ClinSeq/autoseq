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
        self.dbSNP = None
        self.interval_list = None
        self.java_options = '--java-options "-Djava.io.tmpdir={} -Xmx10g"'
        self.output = None
        self.jobname = "gatk-haplotype-germline"

    def command(self):
        haplotypecaller_cmd = "gatk {}".format(self.java_options.format(self.scratch)) +\
                        " HaplotypeCaller " + \
                        required(" -R ", self.reference_sequence) + \
                        required(" -I ", self.input_bam) + \
                        " -L " + self.interval_list + \
                        " --dbsnp " + self.dbSNP + \
                        required(" -O ", self.output) 

        return haplotypecaller_cmd

class StrelkaGermline(Job):
    def __init__(self, input_bam=None, normal_id=None, reference_sequence=None,
                 target_bed=None, output_dir=None, output_filtered_vcf=None ):
        Job.__init__(self)
        self.input_bam = input_bam
        self.normal_id = normal_id
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output_dir = output_dir
        self.output_filtered_vcf = output_filtered_vcf
        
    def command(self):
        required("", self.input_bam)
        required("", self.reference_sequence)

        # configuration
        configure_strelkagermline = "configureStrelkaGermlineWorkflow.py " + \
                                    " --bam " + self.input_bam + \
                                    " --ref " +  self.reference_sequence + \
                                    " --targeted --callRegions " + self.target_bed + \
                                    " --runDir " + self.output_dir
        cmd = configure_strelkagermline + " && " + self.output_dir+"/runWorkflow.py -m local -j 20"

        filter_passed_variants = "zcat " + self.output_dir + "/results/variants/variants.vcf.gz" + \
                                " | awk 'BEGIN { OFS = \"\\t\"} /^#/ { print $0 } {if($7==\"PASS\") print $0 }' " + \
                                " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output_filtered_vcf)
        
        return " && ".join([cmd, filter_passed_variants])

class VarDict(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 reference_dict=None, target_bed=None, output=None, min_alt_frac=0.1, min_num_reads=None,
                 blacklist_bed=None, threads = 1):
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
        self.threads = threads

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        
        mqThres="20"
        bqThres="30"

        freq_filter = (" bcftools filter -e 'STATUS !~ \".*Somatic\" || P0.01Likely=1 || InDelLikely=1' 2> /dev/null "
                       "| %s -c 'from autoseq.util.bcbio import depth_freq_filter_input_stream; import sys; print depth_freq_filter_input_stream(sys.stdin, %s, \"%s\")' " %
                       (sys.executable, 0, 'bwa'))

        somatic_filter = (" sed 's/\\.*Somatic\\\"/Somatic/' "  # changes \".*Somatic\" to Somatic
                          "| sed 's/REJECT,Description=\".*\">/REJECT,Description=\"Not Somatic via VarDict\">/' ")

        blacklist_filter = " | intersectBed -a . -b {} | ".format(self.blacklist_bed)

        cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
              optional("-f ", self.min_alt_frac) + \
              required("-N ", self.tumorid) + \
              optional("-r ", self.min_num_reads) + \
              " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
              optional("-mfreq ", self.min_alt_frac) + \
              optional("-nmfreq ", self.min_alt_frac) + \
              optional("-th ", self.threads) + \
              " -c 1 -S 2 -E 3 -g 4 -Q 10 -B 100 -w 170 -O {} -q {} ".format(mqThres, bqThres) + \
              required("", self.target_bed) + \
              " | testsomatic.R " + \
              " | var2vcf_paired.pl -P 0.05 -m 4.25 -M " + required("-f ", self.min_alt_frac) + \
              " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
              " -q {} -Q {} -d 50 ".format(bqThres, mqThres) + \
              optional("-v ", self.min_num_reads) + \
              " | " + freq_filter + " | " + somatic_filter + " | " + fix_ambiguous_cl() + " | " + remove_dup_cl() + \
              " | vcfstreamsort -w 1000 " + \
              " | " + vt_split_and_leftaln(self.reference_sequence, allow_ref_mismatches=True) + \
              " | bcftools view --apply-filters .,PASS " + \
              " | vcfsorter.pl {} /dev/stdin ".format(self.reference_dict) + \
              conditional(blacklist_filter, self.blacklist_bed) + \
              " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output)
        return cmd

class StrelkaSomatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumor_id=None, normal_id=None, reference_sequence=None,
                 target_bed=None, output_dir=None, output_snvs_vcf=None, output_indels_vcf=None ):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumor_id = tumor_id
        self.normal_id = normal_id
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output_dir = output_dir
        self.output_snvs_vcf = output_snvs_vcf
        self.output_indels_vcf = output_indels_vcf
        
    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration
        configure_strelkasomatic = "configureStrelkaSomaticWorkflow.py --targeted " + \
                                    " --normalBam " + self.input_normal + \
                                    " --tumorBam " + self.input_tumor + \
                                    " --ref " +  self.reference_sequence + \
                                    " --callRegions " + self.target_bed + \
                                    " --runDir " + self.output_dir 
        
        cmd = configure_strelkasomatic + " && " + self.output_dir+"/runWorkflow.py -m local -j 20"

        filter_pass_snvs = "zcat " + self.output_dir + "/results/variants/somatic.snvs.vcf.gz" + \
                      " | awk 'BEGIN { OFS = \"\\t\"} /^#/ { print $0 } {if($7==\"PASS\") print $0 }' " + \
                      " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output_snvs_vcf)

        filter_pass_indels = "zcat " + self.output_dir + "/results/variants/somatic.indels.vcf.gz" + \
                      " | awk 'BEGIN { OFS = \"\\t\"} /^#/ { print $0 } {if($7==\"PASS\") print $0 }' " + \
                      " | bgzip > {output} && tabix -p vcf {output}".format(output=self.output_indels_vcf)

        return " && ".join([cmd, filter_pass_snvs, filter_pass_indels])

class Mutect2Somatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumor_id=None, normal_id=None, reference_sequence=None,
                 target_bed=None, output=None, bamout=None, interval_list=None, output_filtered=None, scratch="/tmp"):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumor_id = tumor_id
        self.normal_id = normal_id
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.output = output
        self.bamout = bamout
        self.interval_list = interval_list
        self.output_filtered = output_filtered
        self.scratch = scratch

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        mqThres="20"
        bqThres="30"

        mutectsomatic_cmd = "gatk --java-options '-Xmx10g -Djava.io.tmpdir=" + self.scratch + "'" + \
                                    " Mutect2 " + \
                                    " -R " +  self.reference_sequence + \
                                    " -I " + self.input_tumor + \
                                    " -I " + self.input_normal + \
                                    " -tumor " + self.tumor_id + \
                                    " -normal " + self.normal_id + \
                                    " -L " + self.interval_list + \
                                    " --base-quality-score-threshold " + bqThres + \
                                    " --callable-depth 50 " + \
                                    " --min-base-quality-score " + bqThres + " --minimum-mapping-quality " + mqThres + \
                                    " --f1r2-median-mq " + mqThres + " --f1r2-min-bq:Integer " + bqThres + \
                                    " -bamout " + self.bamout + \
                                    " -O " + self.output

        filter_mutect_calls = "gatk --java-options '-Xmx10g -Djava.io.tmpdir=" + self.scratch + "'" + \
                                " FilterMutectCalls " + \
                                " -I " + self.bamout + \
                                " -R " +  self.reference_sequence + \
                                " --max-alt-allele-count 2 " + \
                                " --min-median-base-quality " + bqThres + \
                                " --min-median-mapping-quality " + mqThres + \
                                " -V " + self.output + \
                                " -O "  + self.output_filtered

      
        return " && ".join([mutectsomatic_cmd, filter_mutect_calls])

class Varscan2Somatic(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 target_bed=None, normal_pileup=None, tumor_pileup=None, outdir=None,
                 output_somatic_snv=None, output_somatic_indel=None, scratch="/tmp", min_alt_frac=0.1,
                 min_num_reads=None):
        Job.__init__(self)
        self.input_tumor = input_tumor
        self.input_normal = input_normal
        self.tumorid = tumorid
        self.normalid = normalid
        self.reference_sequence = reference_sequence
        self.target_bed = target_bed
        self.normal_pileup = normal_pileup
        self.tumor_pileup = tumor_pileup
        self.outdir = outdir
        self.output_somatic_snv = output_somatic_snv
        self.output_somatic_indel = output_somatic_indel
        self.scratch = scratch
        self.min_alt_frac = min_alt_frac
        self.min_num_reads = min_num_reads

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)
        required("", self.reference_sequence)

        # configuration
        mqThres="20"
        bqThres="30"
        
        normal_mpileup_cmd = "samtools mpileup --no-BAQ --max-depth 50000 --min-MQ " + mqThres + " --min-BQ " + bqThres + " -C50 -f " + \
                             self.reference_sequence + " -l " + self.target_bed + " " + self.input_normal + " > " + \
                             self.normal_pileup
        tumor_mpileup_cmd = "samtools mpileup --no-BAQ --max-depth 50000 --min-MQ " + mqThres + " --min-BQ " + bqThres + " -C50 -f " + \
                            self.reference_sequence + " -l " + self.target_bed + " "  + self.input_tumor + " > " + \
                            self.tumor_pileup
        

        # varscan somatic
        varscan_cmd = "varscan -Xmx10g -Djava.io.tmpdir=" + self.scratch + \
                      " somatic " + self.normal_pileup + " " + self.tumor_pileup + \
                      " --output-snp " + self.outdir + "/varscan.snp.vcf " + \
                      " --output-indel " + self.outdir + "/varscan.indel.vcf " + \
                      optional("--min-var-freq ", self.min_alt_frac) + \
                      " --min-coverage 50 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 0 --output-vcf 1"

        # varscan processSomatic
        process_somatic = "varscan -Xmx10g -Djava.io.tmpdir="+ self.scratch + \
                         " processSomatic " + self.outdir + "/varscan.indel.vcf " + \
                         optional("--min-tumor-freq ", self.min_alt_frac) + \
                         " && varscan -Xmx10g -Djava.io.tmpdir=" + self.scratch + \
                         " processSomatic " + self.outdir + "/varscan.snp.vcf " + \
                         optional("--min-tumor-freq ", self.min_alt_frac)
        
        # create variant list on format that bam-readcount accepts
        variant_lists = "awk 'BEGIN {OFS=\"\\t\"} {if (!/^#/) { isDel=(length($4) > length($5)) ? 1 : 0; print $1,($2+isDel),($2+isDel); }}' " + \
                        self.outdir + "/varscan.indel.Somatic.hc.vcf > " + \
                        self.outdir + "/varscan.indel.Somatic.hc.var && " + \
                        "awk 'BEGIN {OFS=\"\\t\"} {if (!/^#/) { isDel=(length($4) > length($5)) ? 1 : 0; print $1,($2+isDel),($2+isDel); }}' " + \
                        self.outdir + "/varscan.snp.Somatic.hc.vcf > " + \
                        self.outdir + "/varscan.snp.Somatic.hc.var"
        
        # link bam index file to format that bam-readcount accepts
        bai_file_link = "ln -svf " + self.input_tumor.replace(".bam", ".bai") + self.input_tumor + ".bai"
        
        # bam-readcount
        bam_readcount = "bam-readcount -q " + mqThres + " -b " + bqThres + " -l " + \
                        self.outdir + "/varscan.indel.Somatic.hc.var -f " + self.reference_sequence + \
                        " -w 1 " + self.input_tumor + " > " + self.outdir + "/varscan.indel.Somatic.hc.readcount && " + \
                        "bam-readcount -q " + mqThres + " -b " + bqThres + " -l " + \
                        self.outdir + "/varscan.snp.Somatic.hc.var -f " + self.reference_sequence + \
                        " -w 1 " + self.input_tumor + " > " + self.outdir + "/varscan.snp.Somatic.hc.readcount && "
        
        # fpfilter
        fpfilter_indel = "varscan -Xmx10g -Djava.io.tmpdir=" + self.scratch + " fpfilter " + \
                         self.outdir + "/varscan.indel.Somatic.hc.vcf " + \
                         self.outdir + "/varscan.indel.Somatic.hc.readcount " + \
                         " --output-file " + self.output_somatic_indel + \
                         " --filtered-file " + self.output_somatic_indel.replace("filtered.vcf", "failed.vcf") + \
                         optional("--min-var-freq ", self.min_alt_frac) + \
                         optional("--min-var-count ", self.min_num_reads) + \
                         " --min-strandedness 0.1 --min-strand-reads 10 --max-basequal-diff 90 " + \
                         " --min-ref-avgrl 70 --min-var-avgrl 70 --max-ref-mmqs 100 --max-var-mmqs 90 --max-mmqs-diff 90"
        
        fpfilter_snv = "varscan -Xmx10g -Djava.io.tmpdir=" + self.scratch + " fpfilter " + \
                         self.outdir + "/varscan.snp.Somatic.hc.vcf " + \
                         self.outdir + "/varscan.snp.Somatic.hc.readcount " + \
                         " --output-file " + self.output_somatic_snv + \
                         " --filtered-file " + self.output_somatic_snv.replace("filtered.vcf", "failed.vcf") + \
                         optional("--min-var-freq ", self.min_alt_frac) + \
                         optional("--min-var-count ", self.min_num_reads) + \
                         " --min-strandedness 0.1 --min-strand-reads 10 --max-basequal-diff 50 " + \
                         " --min-ref-avgrl 70 --min-var-avgrl 70 --max-ref-mmqs 100 --max-var-mmqs 180 --max-mmqs-diff 180"
        
        # pileup deletion
        delete_pileup = "rm " + self.normal_pileup + " && rm " + self.tumor_pileup
        
        return " && ".join([normal_mpileup_cmd, tumor_mpileup_cmd, varscan_cmd, process_somatic, variant_lists, bai_file_link, bam_readcount, fpfilter_indel, fpfilter_snv, delete_pileup])

class VarDictForPureCN(Job):
    def __init__(self, input_tumor=None, input_normal=None, tumorid=None, normalid=None, reference_sequence=None,
                 reference_dict=None, target_bed=None, output=None, min_alt_frac=0.1, min_num_reads=None, dbsnp=None, threads = 1):
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
        self.threads = threads

    def command(self):
        required("", self.input_tumor)
        required("", self.input_normal)

        mqThres="20"
        bqThres="30"

        tmp_vcf = "{scratch}/{uuid}.vcf.gz".format(scratch=self.scratch, uuid=uuid.uuid4())

        # run vardict without removing non-somatic variants, and adding "SOMATIC" INFO field for somatic variants
        vardict_cmd = "vardict-java " + required("-G ", self.reference_sequence) + \
                      optional("-f ", self.min_alt_frac) + \
                      required("-N ", self.tumorid) + \
                      optional("-r ", self.min_num_reads) + \
                      " -b \"{}|{}\" ".format(self.input_tumor, self.input_normal) + \
                      optional("-mfreq ", self.min_alt_frac) + \
                      optional("-nmfreq ", self.min_alt_frac) + \
                      optional("-th ", self.threads) + \
                      " -c 1 -S 2 -E 3 -g 4 -Q 10 -B 100 -w 170 -O {} -q {} ".format(mqThres, bqThres) + \
                      required("", self.target_bed) + \
                      " | testsomatic.R " + \
                      " | var2vcf_paired.pl -P 0.05 -m 4.25 " + required("-f ", self.min_alt_frac) + \
                      " -N \"{}|{}\" ".format(self.tumorid, self.normalid) + \
                      " -q {} -Q {} -d 50 ".format(bqThres, mqThres) + \
                      optional("-v ", self.min_num_reads) + \
                      " | " + fix_ambiguous_cl() + " | " + remove_dup_cl() + \
                      " | sed 's/Somatic;/Somatic;SOMATIC;/g' " + \
                      " | sed '/^#CHROM/i ##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic event\">' " + \
                      " | vcfstreamsort -w 1000 " + \
                      " | bcftools filter -e 'P0.01Likely=1 || InDelLikely=1' 2> /dev/null " + \
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

class SomaticSeq(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_normal = None
    self.input_tumor = None
    self.reference_sequence = None
    self.input_mutect_vcf = None
    self.input_varscan_snv = None
    self.input_varscan_indel = None
    self.input_vardict_vcf = None
    self.input_strelka_snv = None
    self.input_strelka_indel = None
    self.out_dir = None
    self.out_snv = None
    self.out_indel = None
    self.output_vcf = None
    self.jobname = 'somaticseq-vcf-merging'

  def command(self):

    somatic_seq_env = "source activate somaticseqenv"

    somatic_seq = "run_somaticseq.py " + \
                  " --output-directory " + self.out_dir + \
                  " --genome-reference " + self.reference_sequence +  \
                  " paired " + \
                  " --tumor-bam-file " + self.input_tumor  + \
                  " --normal-bam-file " + self.input_normal + \
                  " --mutect2-vcf " + self.input_mutect_vcf + \
                  " --varscan-snv " + self.input_varscan_snv + \
                  " --varscan-indel " + self.input_varscan_indel + \
                  " --vardict-vcf " + self.input_vardict_vcf + \
                  " --strelka-snv " + self.input_strelka_snv + \
                  " --strelka-indel " + self.input_strelka_indel

    deactivate_ssenv = "source deactivate"

    merge_vcf = "java " + "-Djava.io.tmpdir=" + self.scratch + \
                " -jar /nfs/PROBIO/autoseq-scripts/GenomeAnalysisTK-3.5.jar " + \
                " -T CombineVariants " + \
                " -R " + self.reference_sequence + \
                " --variant " + self.out_snv + \
                " --variant " + self.out_indel + \
                " --assumeIdenticalSamples " + \
                " | bgzip > " + self.output_vcf 
    
    tabix_vcf = "tabix -p vcf {} ".format(self.output_vcf)

    return " && ".join([somatic_seq_env, somatic_seq, deactivate_ssenv, merge_vcf, tabix_vcf])

class VEP(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_vcf = None
        self.reference_sequence = None
        self.vep_dir = None
        self.brca_exchange_vcf = None
        self.jobname = "vep"
        self.additional_options = ""

    def command(self):
        
        fork = ""
        if self.threads > 1:  # vep does not accept "--fork 1", so need to check.
            fork = " --fork {} ".format(self.threads)

        cmdstr = "vep --vcf --output_file STDOUT " + \
                 self.additional_options + required("--dir ", self.vep_dir) + \
                 required("--fasta ", self.reference_sequence) + \
                 required("-i ", self.input_vcf) + \
                 " --check_existing  --total_length --allele_number " + \
                 " --no_escape --no_stats --everything --offline " + \
                 " --custom {},,vcf,exact,0,ClinicalSignificance ".format(self.brca_exchange_vcf) + \
                 fork + " > " + required("", self.output_vcf) 
                 # " && tabix -p vcf {}".format(self.output_vcf)

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

class MergeVCF(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_vcf_hc = None
    self.input_vcf_strelka = None
    self.output_vcf = None
    self.reference_genome = None

  def command(self):

    merge_vcf = "java " + "-Djava.io.tmpdir=" + self.scratch +\
                " -jar /nfs/PROBIO/autoseq-scripts/GenomeAnalysisTK-3.5.jar " + \
                " -T CombineVariants " + \
                " -R " + self.reference_genome + \
                " --variant:haplotypecaller " + self.input_vcf_hc + \
                " --variant:strelka " + self.input_vcf_strelka + \
                " -genotypeMergeOptions PRIORITIZE " + \
                " -priority haplotypecaller,strelka " + \
                " | bgzip > {} ".format(self.output_vcf)
    
    tabix_vcf = "tabix -p vcf {} ".format(self.output_vcf)

    return " && ".join([merge_vcf, tabix_vcf])

class GenerateIGVNavInput(Job):
  def __init__(self):
    Job.__init__(self)
    self.input_vcf = None
    self.oncokb_db = None
    self.output = None
    self.vcftype = None

  def command(self):
    """
    Script to generate tab limited file from final vep annotated vcf file.(both somatic and germline)
    This will annotate the oncogenecity also to final txt file

    usage: generateIGVnavInput.py input.vcf oncokb_allVariants.txt vcftype(somatic or germline) --output ouput.vcf
    """
    return "generateIGVnavInput.py {} {} {} --output {}".format(self.input_vcf, self.oncokb_db, self.vcftype, self.output)

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
        return "vep_install --SPECIES homo_sapiens_vep --AUTO c --ASSEMBLY GRCh37 --NO_HTSLIB " + \
               required("--CACHEDIR ", self.output_dir) + \
               " && vep_convert_cache " + required("--dir ", self.output_dir) + \
               " --species homo_sapiens --version 95_GRCh37"

def call_somatic_variants(pipeline, cancer_bam, normal_bam, cancer_capture, normal_capture,
                          target_name, outdir, callers=['vardict','strelka','mutect2', 'varscan'],
                          min_alt_frac=0.1, min_num_reads=None, scratch="/tmp"):
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
                          output="{}/variants/vardict/{}-{}.vardict-somatic.vcf.gz".format(outdir, cancer_capture_str, normal_capture_str),
                          min_alt_frac=min_alt_frac, min_num_reads=min_num_reads,
                          blacklist_bed=blacklist_bed,
                          threads=pipeline.maxcores
                          )

        vardict.jobname = "vardict/{}".format(cancer_capture_str)
        pipeline.add(vardict)
        d['vardict'] = vardict.output


    if 'strelka' in callers:
        strelka_somatic = StrelkaSomatic(input_tumor=cancer_bam, input_normal=normal_bam, tumor_id=tumor_sample_str,
                          normal_id=normal_sample_str,
                          reference_sequence=pipeline.refdata['reference_genome'],
                          target_bed=pipeline.refdata['targets'][capture_name]['targets-bed-slopped20-gz'],
                          output_dir="{}/variants/{}-{}-strelka-somatic".format(outdir, normal_capture_str, cancer_capture_str),
                          output_snvs_vcf= "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.snvs.vcf.gz".format(outdir, normal_capture_str, cancer_capture_str),
                          output_indels_vcf= "{}/variants/{}-{}-strelka-somatic/results/variants/somatic.passed.indels.vcf.gz".format(outdir, normal_capture_str, cancer_capture_str),
                          )
        strelka_somatic.jobname = "strelka-somatic-workflow/{}".format(cancer_capture_str)
        pipeline.add(strelka_somatic)
        d['strelka_snvs'] = strelka_somatic.output_snvs_vcf
        d['strelka_indels'] = strelka_somatic.output_indels_vcf

    if 'mutect2' in callers:
        mutect_somatic = Mutect2Somatic(input_tumor=cancer_bam, input_normal=normal_bam, tumor_id=tumor_sample_str,
                          normal_id=normal_sample_str,
                          reference_sequence=pipeline.refdata['reference_genome'],
                          output="{}/variants/mutect/{}-{}-gatk-mutect-somatic.vcf.gz".format(outdir, normal_capture_str, cancer_capture_str),
                          bamout="{}/variants/mutect/{}-{}-mutect.bam".format(outdir, normal_capture_str, cancer_capture_str),
                          interval_list=pipeline.refdata['targets'][capture_name]['targets-interval_list-slopped20'],
                          output_filtered="{}/variants/mutect/{}-{}-gatk-mutect-somatic-filtered.vcf.gz".format(outdir, normal_capture_str, cancer_capture_str),
                          scratch = scratch
                          )
        mutect_somatic.jobname = "mutect2-somatic/{}".format(cancer_capture_str)
        pipeline.add(mutect_somatic)
        d['mutect2'] = mutect_somatic.output_filtered

    if 'varscan' in callers:
        varscan_somatic = Varscan2Somatic(input_tumor=cancer_bam, input_normal=normal_bam, tumorid=tumor_sample_str,
                            normalid=normal_sample_str,
                            reference_sequence=pipeline.refdata['reference_genome'],
                            normal_pileup="{}/variants/varscan/{}.pileup".format(outdir, normal_capture_str),
                            tumor_pileup="{}/variants/varscan/{}.pileup".format(outdir, cancer_capture_str),
                            outdir = outdir,
                            target_bed=pipeline.refdata['targets'][capture_name]['targets-bed-slopped20-gz'],
                            output_somatic_snv="{}/variants/varscan/{}-{}-varscan.snp.Somatic.hc.filtered.vcf".format(outdir, normal_capture_str, cancer_capture_str),
                            output_somatic_indel="{}/variants/varscan/{}-{}-varscan.indel.Somatic.hc.filtered.vcf".format(outdir, normal_capture_str, cancer_capture_str),
                            min_alt_frac=min_alt_frac,
                            min_num_reads=min_num_reads,
                            scratch = pipeline.scratch
                            )
        varscan_somatic.jobname = "varscan-somatic/{}".format(cancer_capture_str)
        pipeline.add(varscan_somatic)
        d['varscan_snv'] = varscan_somatic.output_somatic_snv
        d['varscan_indel'] = varscan_somatic.output_somatic_indel

    return d

from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from autoseq.util.path import normpath, stripsuffix
from autoseq.tools.alignment import align_library
from autoseq.tools.cnvcalling import QDNASeq
from autoseq.util.library import find_fastqs
from autoseq.tools.picard import PicardCollectInsertSizeMetrics, PicardCollectOxoGMetrics, \
    PicardMergeSamFiles, PicardMarkDuplicates, PicardCollectHsMetrics, PicardCollectWgsMetrics
from autoseq.tools.variantcalling import Freebayes, VEP, VcfAddSample, call_somatic_variants
from autoseq.tools.intervals import MsiSensor
from autoseq.tools.cnvcalling import CNVkit
from autoseq.tools.contamination import ContEst, ContEstToContamCaveat, CreateContestVCFs
from autoseq.tools.qc import *
from autoseq.util.clinseq_barcode import *
import collections, logging


class SinglePanelResults(object):
    """
    Represents the results generated by performing analysis on a unique sample library capture,
    irrespective of the sample type.
    """
    def __init__(self):
        self.merged_bamfile = None

        # CNV kit outputs:
        self.cnr = None
        self.cns = None
        
        # Coverage QC call:
        self.cov_qc_call = None


class CancerVsNormalPanelResults(object):
    """
    Represents the results generated by performing a paired analysis comparing a cancer and a normal capture.
    """
    def __init__(self):
        self.somatic_vcf = None,
        self.msi_output = None,
        self.hzconcordance_output = None,
        self.vcf_addsample_output = None,
        self.normal_contest_output = None,
        self.cancer_contest_output = None,
        self.cancer_contam_call = None


class ClinseqPipeline(PypedreamPipeline):
    """
    A pipeline for processing clinseq cancer genomics.
    """
    def __init__(self, sampledata, refdata, outdir, libdir, maxcores=1,
                 scratch="/scratch/tmp/tmp", **kwargs):
        """
        :param sampledata: A dictionary specifying the clinseq barcodes of samples of different types.
        :param refdata: A dictionary specifying the reference data used for configuring the pipeline jobs.
        :param outdir: Output folder location string.
        :param libdir: String specifying location of the library fastq files.
        :param maxcores: Maximum number of cores to use concurrently in this analysis.
        :param scratch: String indicating folder in which jobs should output all temporary files.
        :param kwargs: Additional key-word arguments.
        """
        PypedreamPipeline.__init__(self, normpath(outdir), **kwargs)
        self.sampledata = sampledata
        self.refdata = refdata
        self.maxcores = maxcores
        self.libdir = libdir
        self.qc_files = []
        self.scratch = scratch

        # Dictionary linking unique captures to corresponding generic single panel
        # analysis results (SinglePanelResults objects as values):
        self.capture_to_results = collections.defaultdict(SinglePanelResults)

        # Dictionary linking unique normal library capture items to their corresponding
        # germline VCF filenames:
        self.normal_capture_to_vcf = {}

        # Dictionary linking (normal capture, cancer capture) pairings to corresponding
        # cancer library capture analysis results (CancerPanelResults objects as values):
        self.normal_cancer_pair_to_results = collections.defaultdict(CancerVsNormalPanelResults)

    def set_germline_vcf(self, normal_capture, vcf_filename):
        """
        Registers the specified vcf filename for the specified normal capture item,
        for this analysis.

        :param normal_capture: Normal panel capture identifier.
        :param vcf_filename: VCF filename to store.
        """

        self.normal_capture_to_vcf[normal_capture] = vcf_filename

    def get_germline_vcf(self, normal_capture):
        """
        Obtain the germline VCF for the given normal sample capture item.

        :param normal_capture: Named tuple indicating a unique library capture.
        :return: The germline VCF filename for the specified normal capture item, or None
        if this has not been configured.
        """
        if normal_capture in self.normal_capture_to_vcf:
            return self.normal_capture_to_vcf[normal_capture]
        else:
            return None

    def set_capture_bam(self, unique_capture, bam):
        """
        Set the bam file corresponding to the specified unique_capture in this analysis.

        :param unique_capture: A UniqueCapture item. 
        :param bam: The bam filename.
        """
        self.capture_to_results[unique_capture].merged_bamfile = bam

    def set_capture_cnr(self, unique_capture, cnr):
        """
        Record the CNR copy number information (CNV kit output) for the given library capture.

        :param unique_capture: Named tuple indicating unique library capture.
        :param cnr: CNR output filename.
        """
        self.capture_to_results[unique_capture].cnr = cnr 

    def set_capture_cns(self, unique_capture, cns):
        """
        Record the CNS copy number information (CNV kit output) for the given library capture.

        :param unique_capture: Named tuple indicating unique library capture.
        :param cnr: CNS output filename.
        """
        self.capture_to_results[unique_capture].cns = cns

    def get_capture_bam(self, unique_capture):
        """
        Retrieve the bam file corresponding to the specified unique_capture in this analysis.

        :param unique_capture: Named tuple indicating unique library capture.
        :return: The corresponding bam filename, or None if it has not been configured.
        """

        if unique_capture in self.get_all_unique_captures():
            return self.capture_to_results[unique_capture].merged_bamfile
        else:
            return None

    def check_sampledata(self):
        """
        Check this pipeline for validity of the sample data. In particular, check that
        each clinseq barcode has a corresponding fastq file, and if not, then modify
        the pipeline's sampledata by removing that clinseq barcode from the analysis.
        """
        def check_clinseq_barcode_for_data(lib):
            if lib:
                filedir = os.path.join(self.libdir, lib)
                if not os.path.exists(filedir):
                    logging.warn("Dir {} does not exists for {}. Not using library.".format(filedir, lib))
                    return None
                if find_fastqs(lib, self.libdir) == (None, None):
                    logging.warn("No fastq files found for {} in dir {}".format(lib, filedir))
                    return None
            logging.debug("Library {} has data. Using it.".format(lib))
            return lib

        for datatype in ['panel', 'wgs']:
            for sample_type in ['N', 'T', 'CFDNA']:
                clinseq_barcodes_with_data = []
                for clinseq_barcode in self.sampledata[datatype][sample_type]:
                    barcode_checked = check_clinseq_barcode_for_data(clinseq_barcode)
                    if barcode_checked:
                        clinseq_barcodes_with_data.append(barcode_checked)

                self.sampledata[datatype][sample_type] = clinseq_barcodes_with_data

    def get_vep(self):
        """
        Indicates whether the VEP folder has been set for this analysis.

        :return: Boolean.
        """
        vep = False
        if self.refdata['vep_dir']:
            vep = True

        return vep

    def get_all_unique_captures(self):
        """
        Obtain all unique sample library captures in this pipeline instance.

        :return: List of unique capture named tuples. 
        """

        return self.capture_to_results.keys()

    def get_unique_normal_captures(self):
        """
        Obtain tuples for all unique normal sample library captures
        in this pipeline instance - not including "WGS" (no) capture items.

        :return: List of named tuples.
        """

        all_unique_captures = self.get_all_unique_captures()
        return filter(lambda unique_capture: unique_capture.sample_type == "N" and
                                             unique_capture.capture_kit_id != "WG",
                      all_unique_captures)

    def get_unique_wgs(self):
        """
        Obtain all unique cancer sample library WGS items in this pipeline
        instance.

        :return: List of named tuples.
        """

        all_unique_captures = self.get_all_unique_captures()
        return filter(lambda unique_capture: unique_capture.capture_kit_id == "WG",
                      all_unique_captures)

    def get_unique_cancer_captures(self):
        """
        Obtain all unique cancer sample library captures items in this
        pipeline instance - not including "WGS" (no) capture items.

        :return: List of named tuples.
        """

        all_unique_captures = self.get_all_unique_captures()
        return filter(lambda unique_capture: unique_capture.sample_type != "N" and
                                             unique_capture.capture_kit_id != "WG",
                      all_unique_captures)

    def get_prep_kit_name(self, prep_kit_code):
        """
        Convert a two-letter library kit code to the corresponding library kit name.

        :param prep_kit_code: Two-letter library prep code. 
        :return: The library prep kit name.
        """

        # FIXME: Move this information to a config JSON file.
        prep_kit_lookup = {"BN": "BIOO_NEXTFLEX",
                           "KH": "KAPA_HYPERPREP",
                           "TD": "THRUPLEX_DNASEQ",
                           "TP": "THRUPLEX_PLASMASEQ",
                           "TF": "THRUPLEX_FD",
                           "TS": "TRUSEQ_RNA",
                           "NN": "NEBNEXT_RNA",
                           "VI": "VILO_RNA"}

        return prep_kit_lookup[prep_kit_code]

    def get_capture_name(self, capture_kit_code):
        """
        Convert a two-letter capture kit code to the corresponding capture kit name.

        :param capture_kit_code: The two-letter capture kit code.
        :return: The capture-kit name.
        """

        # FIXME: Move this information to a config JSON file.
        capture_kit_loopkup = {"CS": "clinseq_v3_targets",
                               "CZ": "clinseq_v4",
                               "EX": "EXOMEV3",
                               "EO": "EXOMEV1",
                               "RF": "fusion_v1",
                               "CC": "core_design",
                               "CD": "discovery_coho",
                               "CB": "big_design",
                               "AL": "alascca_targets",
                               "TT": "test-regions",
                               "CP": "progression",
                               "CM": "monitor"
                               }

        if capture_kit_code == 'WG':
            return 'lowpass_wgs'

        else:
            return capture_kit_loopkup[capture_kit_code]

    def get_all_clinseq_barcodes(self):
        """
        :return: All clinseq barcodes included in this clinseq analysis pipeline's panel data.
        """
        
        all_clinseq_barcodes = \
            self.sampledata['T'] + \
            self.sampledata['N'] + \
            self.sampledata['CFDNA']
        return filter(lambda bc: bc != None, all_clinseq_barcodes)

    def get_unique_capture_to_clinseq_barcodes(self):
        """
        Retrieves all clinseq barcodes for this clinseq analysis, and organises them according
        to unique library captures.

        :return: A dictionary with tuples indicating unique library captures as keys,
        and barcode lists as values.
        """

        capture_to_barcodes = collections.defaultdict(list)
        for clinseq_barcode in self.get_all_clinseq_barcodes():
            unique_capture = parse_capture_tuple(clinseq_barcode)
            capture_to_barcodes[unique_capture].append(clinseq_barcode)

        return capture_to_barcodes

    def merge_and_rm_dup(self, unique_capture, input_bams):
        """
        Configures Picard merging and duplicate marking, for the specified group input bams,
        which should all correspond to the specified sample library capture.
        
        Registers the final output bam file for this library capture in this analysis.

        :param sample_type: Clinseq sample type
        :param sample_id: Sample ID
        :param prep_kit_id: Two-letter prep kit ID
        :param capture_kit_id: Two-letter capture kit ID
        :input_bams: The bam filenames for which to do merging and duplicate marking
        """

        # Strings indicating the sample and capture, for use in output file names below:
        sample_str = "{}-{}".format(unique_capture.sample_type, unique_capture.sample_id)
        capture_str = "{}-{}-{}".format(sample_str, unique_capture.prep_kit_id, unique_capture.capture_kit_id)

        # Configure merging:
        merged_bam_filename = \
            "{}/bams/{}/{}.bam".format(self.outdir, unique_capture.capture_kit_id, capture_str)
        merge_bams = PicardMergeSamFiles(input_bams, merged_bam_filename)
        merge_bams.is_intermediate = True
        merge_bams.jobname = "picard-mergebams-{}".format(sample_str)
        self.add(merge_bams)

        # Configure duplicate marking:
        mark_dups_bam_filename = \
            "{}/bams/{}/{}-nodups.bam".format(self.outdir, unique_capture.capture_kit_id, capture_str)
        mark_dups_metrics_filename = \
            "{}/qc/picard/{}/{}-markdups-metrics.txt".format(
                self.outdir, unique_capture.capture_kit_id, capture_str)
        markdups = PicardMarkDuplicates(
            merge_bams.output_bam, mark_dups_bam_filename, mark_dups_metrics_filename)
        markdups.is_intermediate = False
        self.add(markdups)

        self.set_capture_bam(unique_capture, markdups.output_bam)

        self.qc_files.append(markdups.output_metrics)

    def configure_fastq_qcs(self):
        """
        Configure QC on all fastq files that exist for this pipeline instance.

        :return: List of qc output filenames.
        """

        qc_files = []
        for clinseq_barcode in self.get_all_clinseq_barcodes():
            curr_fqs = find_fastqs(clinseq_barcode, self.libdir)
            for fq in curr_fqs:
                fastqc = FastQC()
                fastqc.input = fq
                fastqc.outdir = "{}/qc/fastqc/".format(self.outdir)
                fastqc.output = "{}/qc/fastqc/{}_fastqc.zip".format(
                    self.outdir, clinseq_barcode)
                fastqc.jobname = "fastqc-{}".format(clinseq_barcode)
                qc_files.append(fastqc.output)
                self.add(fastqc)

        return qc_files

    def configure_align_and_merge(self):
        """
        Configure the aligning of the fastq files for all clinseq barcodes in this pipeline,
        and configure merging of the resulting bam files organised according to unique
        sample library captures (including "WGS" captures - i.e. no capture).
        """

        capture_to_barcodes = self.get_unique_capture_to_clinseq_barcodes()
        for unique_capture in capture_to_barcodes.keys():
            curr_bamfiles = []
            capture_kit = unique_capture.capture_kit_id
            for clinseq_barcode in capture_to_barcodes[unique_capture]:
                curr_bamfiles.append(
                    align_library(self,
                                  fq1_files=find_fastqs(clinseq_barcode, self.libdir)[0],
                                  fq2_files=find_fastqs(clinseq_barcode, self.libdir)[1],
                                  clinseq_barcode=clinseq_barcode,
                                  ref=self.refdata['bwaIndex'],
                                  outdir= "{}/bams/{}".format(self.outdir, capture_kit),
                                  maxcores=self.maxcores))

            self.merge_and_rm_dup(unique_capture, curr_bamfiles)

    def call_germline_variants(self, normal_capture, bam):
        """
        Configure calling of germline variants for a normal sample library capture,
        and configure VEP if specified in the analysis.

        :param normal_capture: The normal sample library capture identifier.
        :param bam: Bam filename input to variant calling.
        """

        targets = self.get_capture_name(normal_capture.capture_kit_id)
        capture_str = "{}-{}-{}".format(normal_capture.sample_id,
                                        normal_capture.prep_kit_id,
                                        normal_capture.capture_kit_id)

        freebayes = Freebayes()
        freebayes.input_bams = [bam]
        freebayes.somatic_only = False
        freebayes.params = None
        freebayes.reference_sequence = self.refdata['reference_genome']
        freebayes.target_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        freebayes.threads = self.maxcores
        freebayes.scratch = self.scratch
        freebayes.output = "{}/variants/{}.freebayes-germline.vcf.gz".format(self.outdir, capture_str)
        freebayes.jobname = "freebayes-germline-{}".format(capture_str)
        self.add(freebayes)

        if self.refdata['vep_dir']:
            vep_freebayes = VEP()
            vep_freebayes.input_vcf = freebayes.output
            vep_freebayes.threads = self.maxcores
            vep_freebayes.reference_sequence = self.refdata['reference_genome']
            vep_freebayes.vep_dir = self.refdata['vep_dir']
            vep_freebayes.output_vcf = "{}/variants/{}.freebayes-germline.vep.vcf.gz".format(self.outdir, capture_str)
            vep_freebayes.jobname = "vep-freebayes-germline-{}".format(capture_str)
            self.add(vep_freebayes)

            self.set_germline_vcf(normal_capture, vep_freebayes.output_vcf)
        else:
            self.set_germline_vcf(normal_capture, freebayes.output)

    def configure_panel_analysis_with_normal(self, normal_capture):
        """
        Configure panel analyses focused on a specific unique normal library capture.
        """

        if normal_capture.sample_type is not "N":
            raise ValueError("Invalid input tuple: " + normal_capture)

        normal_bam = self.get_capture_bam(normal_capture)
        # Configure germline variant calling:
        self.call_germline_variants(normal_capture, normal_bam)

        # For each unique cancer library capture, configure a comparative analysis against
        # this normal capture:
        for cancer_capture in self.get_unique_cancer_captures():
            self.configure_panel_analysis_cancer_vs_normal(
                normal_capture, cancer_capture)

    def configure_single_capture_analysis(self, unique_capture):
        """
        Configure all general analyses to perform given a single sample library capture.
        """

        input_bam = self.get_capture_bam(unique_capture)
        sample_str = compose_sample_str(unique_capture)
        targets = self.get_capture_name(unique_capture.capture_kit_id)

        # Configure CNV kit analysis:
        cnvkit = CNVkit(input_bam=input_bam,
                        output_cnr="{}/cnv/{}.cnr".format(self.outdir, sample_str),
                        output_cns="{}/cnv/{}.cns".format(self.outdir, sample_str),
                        scratch=self.scratch)

        # If we have a CNVkit reference
        if self.refdata['targets'][targets]['cnvkit-ref']:
            cnvkit.reference = self.refdata['targets'][targets]['cnvkit-ref']
        else:
            cnvkit.targets_bed = self.refdata['targets'][targets]['targets-bed-slopped20']

        cnvkit.jobname = "cnvkit/{}".format(sample_str)

        # Register the result of this analysis:
        self.set_capture_cnr(unique_capture, cnvkit.output_cnr)
        self.set_capture_cns(unique_capture, cnvkit.output_cns)

        self.add(cnvkit)

    def configure_lowpass_analyses(self):
        """
        Configure generic analyses of all low-pass whole-genome sequencing
        data for this clinseq pipeline, under the assumption that alignment and
        bam file merging has already been performed."""

        for unique_wgs in self.get_unique_wgs():
            self.configure_single_wgs_analyses(unique_wgs)

    def configure_single_wgs_analyses(self, unique_wgs):
        """
        Configure generic analyses of a single WGS item in the pipeline.

        :param unique_wgs: An identifier for a single unique library WGS.
        """

        input_bam = self.get_capture_bam(unique_wgs)
        sample_str = compose_sample_str(unique_wgs)

        qdnaseq = QDNASeq(input_bam,
                          output_segments="{}/cnv/{}-qdnaseq.segments.txt".format(
                              self.outdir, sample_str),
                          background=None
                          )
        self.add(qdnaseq)

    def run_wgs_bam_qc(self, bams):
        """
        Run QC on wgs bams
        :param bams: list of bams
        :return: list of generated files
        """
        qc_files = []
        logging.debug("bams are {}".format(bams))
        for bam in bams:
            basefn = stripsuffix(os.path.basename(bam), ".bam")
            isize = PicardCollectInsertSizeMetrics()
            isize.input = bam
            isize.jobname = "picard-isize-{}".format(basefn)
            isize.output_metrics = "{}/qc/picard/wgs/{}.picard-insertsize.txt".format(self.outdir, basefn)
            self.add(isize)

            wgsmetrics = PicardCollectWgsMetrics()
            wgsmetrics.input = bam
            wgsmetrics.reference_sequence = self.refdata['reference_genome']
            wgsmetrics.output_metrics = "{}/qc/picard/wgs/{}.picard-wgsmetrics.txt".format(self.outdir, basefn)
            wgsmetrics.jobname = "picard-wgsmetrics-{}".format(basefn)
            self.add(wgsmetrics)

            qc_files += [isize.output_metrics, wgsmetrics.output_metrics]

        return qc_files

    def configure_panel_analyses(self):
        """
        Configure generic analyses of all panel data for this clinseq pipeline,
        assuming that alignment and bam file merging has been performed.
        """

        # Configure analyses to be run on all unique panel captures individually:
        for unique_capture in self.get_unique_cancer_captures():
            self.configure_single_capture_analysis(unique_capture)

        # Configure a separate group of analyses for each unique normal library capture:
        for normal_capture in self.get_unique_normal_captures():
            self.configure_panel_analysis_with_normal(normal_capture)

    def configure_somatic_calling(self, normal_capture, cancer_capture):
        """
        Configure somatic variant calling in this pipeline, for a specified pairing
        of normal and cancer library capture events. 

        :param normal_capture: Named tuple indicating normal library capture.
        :param cancer_capture: Named tuple indicating cancer library capture.
        """
        # FIXME: Need to fix the configuration of the min_alt_frac threshold, rather than hard-coding it here:
        somatic_variants = call_somatic_variants(
            self, tbam=self.get_capture_bam(normal_capture), nbam=self.get_capture_bam(cancer_capture),
            tlib=compose_sample_str(cancer_capture), nlib=compose_sample_str(normal_capture),
            target_name=self.get_capture_name(cancer_capture.capture_kit_id),
            refdata=self.refdata, outdir=self.outdir,
            callers=['vardict'], vep=self.get_vep(), min_alt_frac=0.02)
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].somatic_vcf = somatic_variants

    def configure_vcf_add_sample(self, normal_capture, cancer_capture):
        """
        Configure VCF updating in this pipeline, for a specified pairing
        of normal and cancer library capture events. 

        :param normal_capture: Named tuple indicating normal library capture.
        :param cancer_capture: Named tuple indicating cancer library capture.
        """

        # Configure VCF add sample:
        vcfaddsample = VcfAddSample()
        vcfaddsample.input_bam = self.get_capture_bam(cancer_capture)
        vcfaddsample.input_vcf = self.get_germline_vcf(normal_capture)
        normal_sample_str = compose_sample_str(normal_capture)
        cancer_sample_str = compose_sample_str(cancer_capture)
        vcfaddsample.samplename = cancer_sample_str
        vcfaddsample.filter_hom = True
        vcfaddsample.output = "{}/variants/{}-and-{}.germline-variants-with-somatic-afs.vcf.gz".format(
            self.outdir, normal_sample_str, cancer_sample_str)
        vcfaddsample.jobname = "vcf-add-sample-{}".format(cancer_sample_str)
        self.add(vcfaddsample)
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vcf_addsample_output = \
            vcfaddsample.output

    def configure_msi_sensor(self, normal_capture, cancer_capture):
        """
        Configure MSI sensor in this pipeline, for a specified pairing
        of normal and cancer library capture events. 

        :param normal_capture: Named tuple indicating normal library capture.
        :param cancer_capture: Named tuple indicating cancer library capture.
        """

        # Configure MSI sensor:
        msisensor = MsiSensor()
        msisensor.msi_sites = self.refdata['targets'][cancer_capture.capture_kit_id]['msisites']
        msisensor.input_normal_bam = self.get_capture_bam(normal_capture)
        msisensor.input_tumor_bam = self.get_capture_bam(cancer_capture)
        normal_capture_str = compose_sample_str(normal_capture)
        cancer_capture_str = compose_sample_str(cancer_capture)
        msisensor.output = "{}/msisensor-{}-{}.tsv".format(
            self.outdir, normal_capture_str, cancer_capture_str)
        msisensor.threads = self.maxcores
        msisensor.jobname = "msisensor-{}-{}".format(normal_capture_str, cancer_capture_str)
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].msi_output = \
            msisensor.output
        self.add(msisensor)

    def configure_hz_conc(self, normal_capture, cancer_capture):
        """
        Configure heterozygote concordance calculation in this pipeline, for a
        specified pairing of normal and cancer library capture events. 

        :param normal_capture: Named tuple indicating normal library capture.
        :param cancer_capture: Named tuple indicating cancer library capture.
        """

        # Configure heterozygote concordance:
        hzconcordance = HeterzygoteConcordance()
        hzconcordance.input_vcf = self.get_germline_vcf(normal_capture)
        hzconcordance.input_bam = self.get_capture_bam(cancer_capture)
        hzconcordance.reference_sequence = self.refdata['reference_genome']
        hzconcordance.target_regions = \
            self.refdata['targets'][cancer_capture.capture_kit_id]['targets-interval_list-slopped20']
        hzconcordance.normalid = compose_sample_str(normal_capture)
        hzconcordance.filter_reads_with_N_cigar = True
        hzconcordance.jobname = "hzconcordance-{}".format(compose_sample_str(cancer_capture))
        hzconcordance.output = "{}/bams/{}-{}-hzconcordance.txt".format(
            self.outdir, compose_sample_str(cancer_capture), compose_sample_str(normal_capture))
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].hzconcordance_output = \
            hzconcordance.output
        self.add(hzconcordance)

    def configure_contest_vcf_generation(self, normal_capture, cancer_capture):
        """
        Configure generation of a contest VCF input file in this pipeline, for a
        specified pairing of normal and cancer library capture events. 

        :param normal_capture: Named tuple indicating normal library capture.
        :param cancer_capture: Named tuple indicating cancer library capture.
        """

        contest_vcf_generation = CreateContestVCFs()
        contest_vcf_generation.input_population_vcf = self.refdata['swegene_common']
        contest_vcf_generation.input_target_regions_bed_1 = normal_capture
        contest_vcf_generation.input_target_regions_bed_2 = cancer_capture
        normal_capture_str = compose_sample_str(normal_capture)
        cancer_capture_str = compose_sample_str(cancer_capture)
        contest_vcf_generation.output = "{}/contamination/pop_vcf_{}-{}.vcf".format(
            normal_capture_str, cancer_capture_str)
        contest_vcf_generation.jobname = "contest_pop_vcf_{}-{}".format(
            normal_capture_str, cancer_capture_str)
        self.add(contest_vcf_generation)
        return contest_vcf_generation.output

    def configure_contest(self, library_capture_1, library_capture_2, contest_vcf):
        """
        Configure running of ContEst in this pipeline, for a specified pairing
        of library capture events. Estimates contamination in the bam file for the
        first library capture, using the bam file for the second library capture as
        a reference comparison.

        :param library_capture_1: Named tuple indicating first library capture.
        :param library_capture_2: Named tuple indicating second library capture.
        :param contest_vcf: Contest population allele frequency VCF input file.
        """

        contest = ContEst()
        contest.reference_genome = self.refdata['reference_genome']
        contest.input_eval_bam = self.get_capture_bam(library_capture_1)
        contest.input_genotype_bam = self.get_capture_bam(library_capture_2)
        contest.input_population_af_vcf = contest_vcf
        # TODO: Is it necessary to create the output subdir contamination somewhere? Check how it's done for e.g. cnvkit.
        contest.output = "{}/contamination/{}.contest.txt".format(self.outdir, compose_sample_str(library_capture_1)) # TODO: Should the analysis id also be in name of out file?
        contest.jobname = "contest_tumor/{}".format(compose_sample_str(library_capture_1))  # TODO: Is it ok that the job name does not contain analysis id, i.e. may not be unique?
        self.add(contest)
        return contest.output

    def configure_contam_qc_call(self, contest_output):
        """
        Configure generation of a contamination QC call in this pipeline,
        based on the specified contest output. Returns the resulting QC output
        filename.

        :param contest_output: ContEst output filename.
        """

        process_contest = ContEstToContamCaveat()
        process_contest.input_contest_results = contest_output
        process_contest.output = "{}/qc/{}-contam-qc-call.json".format(self.outdir, self.sampledata['panel']['T'])
        self.add(process_contest)
        return process_contest.output

    def configure_contamination_estimate(self, normal_capture, cancer_capture):
        """
        Configure contamination estimatates for a given normal, cancer library capture
        pairing.
        
        :param normal_capture: Namedtuple indicating a normal library capture.
        :param cancer_capture: Namedtuple indicating a cancer library capture. 
        """
        # Configure generation of the contest VCF input file:
        intersection_contest_vcf = \
            self.configure_contest_vcf_generation(normal_capture, cancer_capture)

        # Configure contest for calculating contamination in the cancer sample:
        cancer_vs_normal_contest_output = \
            self.configure_contest(cancer_capture, normal_capture, intersection_contest_vcf)

        # Configure contest for calculating contamination in the normal sample:
        normal_vs_cancer_contest_output = \
            self.configure_contest(normal_capture, cancer_capture, intersection_contest_vcf)

        # Configure cancer sample contamination QC call:
        cancer_contam_call = self.configure_contam_qc_call(self, cancer_vs_normal_contest_output)

        # Register the outputs of running contest:
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].normal_contest_output = \
            normal_vs_cancer_contest_output
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].cancer_contest_output = \
            cancer_vs_normal_contest_output
        self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].cancer_contam_call = \
            cancer_contam_call

    def configure_panel_analysis_cancer_vs_normal(self, normal_capture, cancer_capture):
        """
        Configures standard paired cancer vs normal panel analyses for the specified unique
        normal and cancer library captures.

        Comprises the following analyses:
        - Somatic variant calling
        - Updating of the germline VCF to take into consideration the cancer sample
        - MSI sensor
        - Heterozygote concordance of the sample pair
        - Contamination estimate of cancer compared with normal and vice versa

        :param normal_capture: A unique normal sample library capture
        :param cancer_capture: A unique cancer sample library capture
        """

        self.configure_somatic_calling(normal_capture, cancer_capture)
        self.configure_vcf_add_sample(normal_capture, cancer_capture)
        self.configure_msi_sensor(normal_capture, cancer_capture)
        self.configure_hz_conc(normal_capture, cancer_capture)
        self.configure_contamination_estimate(normal_capture, cancer_capture)

    def configure_all_lowpass_qcs(self):
        """
        Configure QC checks for all low-pass whole genome data in this pipeline.
        """

        for unique_wgs in self.get_unique_wgs():
            self.qc_files += \
                self.configure_wgs_qc(unique_wgs)

    def configure_wgs_qc(self, unique_wgs):
        """
        Configure QC checks for the specified unique WGS item in this pipeline.

        :param unique_wgs: A named tuple identifying a single unique WGS item.
        :return: QC files output files resulting from the QC analysis configuration.
        """

        bam = self.get_capture_bam(unique_wgs)

        qc_files = []

        basefn = stripsuffix(os.path.basename(bam), ".bam")
        isize = PicardCollectInsertSizeMetrics()
        isize.input = bam
        isize.jobname = "picard-isize-{}".format(basefn)
        isize.output_metrics = "{}/qc/picard/wgs/{}.picard-insertsize.txt".format(self.outdir, basefn)
        self.add(isize)

        wgsmetrics = PicardCollectWgsMetrics()
        wgsmetrics.input = bam
        wgsmetrics.reference_sequence = self.refdata['reference_genome']
        wgsmetrics.output_metrics = "{}/qc/picard/wgs/{}.picard-wgsmetrics.txt".format(self.outdir, basefn)
        wgsmetrics.jobname = "picard-wgsmetrics-{}".format(basefn)
        self.add(wgsmetrics)

        qc_files += [isize.output_metrics, wgsmetrics.output_metrics]

        return qc_files

    def configure_all_panel_qcs(self):
        """
        Configures QC checks for all panel data in this pipeline.
        """

        for unique_capture in self.get_all_unique_captures():
            self.qc_files += \
                self.configure_panel_qc(unique_capture)

    def configure_multi_qc(self):
        """
        Configures MultiQC for this pipeline. self.qc_files must be fully populated
        in order for MultiQC to use all relevant input files.
        """

        multiqc = MultiQC()
        multiqc.input_files = self.qc_files
        multiqc.search_dir = self.outdir
        multiqc.output = "{}/multiqc/{}-multiqc".format(self.outdir, self.sampledata['sdid'])
        multiqc.jobname = "multiqc-{}".format(self.sampledata['sdid'])
        self.add(multiqc)

    def configure_panel_qc(self, unique_capture):
        """
        Configure QC analyses for a given library capture.

        :param unique_capture: Named tuple identifying a sample library capture.
        :return: list of QC output files for this capture.
        """

        bam = self.get_capture_bam(unique_capture)

        targets = self.get_capture_name(unique_capture.capture_kit_id)
        logging.debug("Adding QC jobs for {}".format(bam))

        capture_str = compose_sample_str(unique_capture)

        isize = PicardCollectInsertSizeMetrics()
        isize.input = bam
        isize.output_metrics = "{}/qc/picard/{}/{}.picard-insertsize.txt".format(
            self.outdir, unique_capture.capture_kit_id, capture_str)
        isize.jobname = "picard-isize-{}".format(capture_str)
        self.add(isize)

        oxog = PicardCollectOxoGMetrics()
        oxog.input = bam
        oxog.reference_sequence = self.refdata['reference_genome']
        oxog.output_metrics = "{}/qc/picard/{}/{}.picard-oxog.txt".format(
            self.outdir, unique_capture.capture_kit_id, capture_str)
        oxog.jobname = "picard-oxog-{}".format(capture_str)
        self.add(oxog)

        hsmetrics = PicardCollectHsMetrics()
        hsmetrics.input = bam
        hsmetrics.reference_sequence = self.refdata['reference_genome']
        hsmetrics.target_regions = self.refdata['targets'][targets][
            'targets-interval_list-slopped20']
        hsmetrics.bait_regions = self.refdata['targets'][targets][
            'targets-interval_list-slopped20']
        hsmetrics.bait_name = targets
        hsmetrics.output_metrics = "{}/qc/picard/{}/{}.picard-hsmetrics.txt".format(
            self.outdir, unique_capture.capture_kit_id, capture_str)
        hsmetrics.jobname = "picard-hsmetrics-{}".format(capture_str)
        self.add(hsmetrics)

        sambamba = SambambaDepth()
        sambamba.targets_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        sambamba.input = bam
        sambamba.output = "{}/qc/sambamba/{}.sambamba-depth-targets.txt".format(
            self.outdir, capture_str)
        sambamba.jobname = "sambamba-depth-{}".format(capture_str)
        self.add(sambamba)

        coverage_hist = CoverageHistogram()
#        if 'alascca_targets' in self.refdata['targets']:
#            alascca_coverage_hist.input_bed = self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
#        else:
        coverage_hist.input_bed = self.refdata['targets'][targets]['targets-bed-slopped20']
        coverage_hist.input_bam = bam
        coverage_hist.output = "{}/qc/{}.coverage-histogram.txt".format(
            self.outdir, capture_str)
        coverage_hist.jobname = "alascca-coverage-hist/{}".format(capture_str)
        self.add(coverage_hist)

        coverage_qc_call = CoverageCaveat()
        coverage_qc_call.input_histogram = coverage_hist.output
        coverage_qc_call.output = "{}/qc/{}.coverage-qc-call.json".format(self.outdir, capture_str)
        coverage_qc_call.jobname = "coverage-qc-call/{}".format(capture_str)
        self.add(coverage_qc_call)
        self.capture_to_results[unique_capture].cov_qc_call = coverage_qc_call.output

        return [isize.output_metrics, oxog.output_metrics, hsmetrics.output_metrics,
                sambamba.output, coverage_hist.output, coverage_qc_call.output]

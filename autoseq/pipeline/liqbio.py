from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.cnvcalling import LiqbioCNAPlot
from autoseq.util.clinseq_barcode import *
from autoseq.tools.structuralvariants import Svcaller, Sveffect, MantaSomaticSV
from autoseq.tools.structuralvariants import SViCT, Svaba, Lumpy, GenerateIGVNavInputSV, Gridss
from autoseq.tools.structuralvariants import AnnotateSvaba
from autoseq.tools.umi import *
from autoseq.tools.alignment import fq_trimming, Realignment
from autoseq.util.library import find_fastqs

__author__ = 'thowhi'


class LiqBioPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, job_params, outdir, libdir, umi, maxcores=1, scratch="/scratch/tmp/tmp",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, job_params, outdir, libdir, umi,
                                 maxcores, scratch, **kwargs)

        # Set the min alt frac value:
        self.default_job_params["vardict-min-alt-frac"] = 0.01
        self.default_job_params["vardict-min-num-reads"] = None
        self.default_job_params["vep-additional-options"] = " --pick --filter_common "

        # Remove clinseq barcodes for which data is not available:
        self.check_sampledata()

        if umi:
            # Configure the umi processes from fastq to bam file:
            self.configure_umi_processing()
        else:
            # Configure alignment and merging of fastq data for all clinseq barcodes:
            self.configure_align_and_merge()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure liqbio-specific panel analyses:
        self.configure_panel_analyses_liqbio(umi)

        # Configure additional msings analysis:
        self.configure_panel_msings_analyses()

        # Configure QC of all panel data:
        self.configure_all_panel_qcs()

        # Configure fastq QCs:
        self.configure_fastq_qcs()

        # Configure the low-pass whole genome analysis:
        self.configure_lowpass_analyses()

        # Configure low-pass whole genome data QC:
        self.configure_all_lowpass_qcs()

        # Configure MultiQC:
        self.configure_multi_qc()
        
        # Configure QC overview plot
        for normal_capture in self.get_mapped_captures_normal():
            for cancer_capture in self.get_mapped_captures_cancer():
                self.configure_qc_overview_plot(normal_capture, cancer_capture)

    def configure_single_capture_analysis_liqbio(self, unique_capture):
        input_bam = self.get_capture_bam(unique_capture, umi=False)
        sample_str = compose_lib_capture_str(unique_capture)

        # Configure svcaller analysis for each event type:
        for event_type in ["DEL", "DUP", "INV", "TRA"]:
            svcaller = Svcaller()
            svcaller.input_bam = input_bam
            svcaller.event_type = event_type
            svcaller.output_bam = "{}/svs/{}-{}.bam".format(self.outdir, sample_str, event_type)
            svcaller.output_gtf = "{}/svs/{}-{}.gtf".format(self.outdir, sample_str, event_type)
            svcaller.reference_sequence = self.refdata["reference_genome"]
            svcaller.scratch = self.scratch
            self.add(svcaller)

            self.set_capture_svs(unique_capture, event_type, (svcaller.output_bam, svcaller.output_gtf))

        # FIXME: This code is kind of nasty, as the self.capture_to_results data structure is
        # getting "pushed too far" in it's usage:
        sveffect = Sveffect()
        sveffect.input_del_gtf = self.capture_to_results[unique_capture].svs["DEL"][1]
        sveffect.input_dup_gtf = self.capture_to_results[unique_capture].svs["DUP"][1]
        sveffect.input_inv_gtf = self.capture_to_results[unique_capture].svs["INV"][1]
        sveffect.input_tra_gtf = self.capture_to_results[unique_capture].svs["TRA"][1]
        sveffect.ts_regions = self.refdata["ts_regions"]
        sveffect.ar_regions = self.refdata["ar_regions"]
        sveffect.fusion_regions = self.refdata["fusion_regions"]
        sveffect.output_combined_bed = "{}/svs/{}_combined.bed".format(self.outdir, sample_str)
        sveffect.output_effects_json = "{}/svs/{}_effects.json".format(self.outdir, sample_str)

        self.add(sveffect)

        self.set_capture_sveffect(unique_capture, sveffect.output_effects_json)

        svcaller_igvinput = GenerateIGVNavInputSV()
        svcaller_igvinput.input_del_gtf = self.capture_to_results[unique_capture].svs["DEL"][1]
        svcaller_igvinput.input_dup_gtf = self.capture_to_results[unique_capture].svs["DUP"][1]
        svcaller_igvinput.input_inv_gtf = self.capture_to_results[unique_capture].svs["INV"][1]
        svcaller_igvinput.input_tra_gtf = self.capture_to_results[unique_capture].svs["TRA"][1]
        svcaller_igvinput.input_vcf = "{}/svs/".format(self.outdir)
        svcaller_igvinput.sdid = sample_str
        svcaller_igvinput.vcftype = 'normal' if '-N-'in sample_str else 'somatic'
        svcaller_igvinput.tool = 'svcaller'
        svcaller_igvinput.output = "{}/svs/igv/".format(self.outdir)
        svcaller_igvinput.jobname = "generate-igvnav-input-svcaller-" + sample_str

        self.add(svcaller_igvinput)


    def configure_panel_analyses_liqbio(self, umi):
        # Configure liqbio analyses to be run on all unique panel captures individually:
        for unique_capture in self.get_mapped_captures_no_wgs():
            self.configure_single_capture_analysis_liqbio(unique_capture)
            self.configure_svict(unique_capture)
            #self.configure_gridss(unique_capture)

        # Configure a liqbio analyses for each normal-cancer pairing:
        for normal_capture in self.get_mapped_captures_normal():
            for cancer_capture in self.get_mapped_captures_cancer():
                self.configure_panel_analysis_cancer_vs_normal_liqbio(
                    normal_capture, cancer_capture, umi)
                self.configure_gridss(normal_capture, cancer_capture)

    def configure_gridss(self, normal_capture, cancer_capture):
        input_normal = self.get_capture_bam(normal_capture, umi=False)
        input_cancer = self.get_capture_bam(cancer_capture, umi=False)

        normal_str = compose_lib_capture_str(normal_capture)
        cancer_str = compose_lib_capture_str(cancer_capture)

        gridss = Gridss()
        gridss.input_normal = input_normal
        gridss.input_tumor = input_cancer
        gridss.reference_sequence = self.refdata['bwaIndex']
        gridss.pondir = self.refdata["pondir"]
        gridss.assembly = "{}/svs/gridss/{}-{}-assembly.bam".format(self.outdir, normal_str, cancer_str)
        gridss.steps = " ALL "
        gridss.workingdir = "{}/svs/gridss/".format(self.outdir)
        gridss.output_full = "{}/svs/gridss/{}-{}-gridss.vcf.gz".format(self.outdir, normal_str, cancer_str)
        gridss.output_filter = "{}/svs/gridss/{}-{}-gridss_filtered.vcf.bgz".format(self.outdir, normal_str, cancer_str)
        gridss.threads = self.maxcores
        gridss.jobname = "gridss-sv-calling-{}-{}".format(normal_str, cancer_str)

        self.add(gridss)

        # gridss_igvinput = GenerateIGVNavInputSV()
        # gridss_igvinput.input_vcf = gridss.output
        # gridss_igvinput.sdid = unique_capture.sdid
        # gridss_igvinput.vcftype = 'normal' if '-N-'in sample_str else 'somatic'
        # gridss_igvinput.tool = 'gridss'
        # gridss_igvinput.output = "{}/svs/igv/{}".format(self.outdir, sample_str)
        # gridss_igvinput.jobname = "generate-igvnav-input-gridss"

        # self.add(gridss_igvinput)

    def configure_svict(self, unique_capture):

        input_bam = self.get_capture_bam(unique_capture, umi=False)
        sample_str = compose_lib_capture_str(unique_capture)

        svict = SViCT()
        svict.input_bam = input_bam
        svict.reference_sequence = self.refdata["reference_genome"]
        svict.output = "{}/svs/svict/{}-svict".format(self.outdir, sample_str)
        svict.output_vcf = "{}/svs/svict/{}-svict.vcf".format(self.outdir, sample_str)
        svict.threads = self.maxcores

        self.add(svict)

        svict_igvinput = GenerateIGVNavInputSV()
        svict_igvinput.input_vcf = svict.output_vcf
        svict_igvinput.sdid = unique_capture.sdid
        svict_igvinput.vcftype = 'normal' if '-N-'in sample_str else 'somatic'
        svict_igvinput.tool = 'svict'
        svict_igvinput.output = "{}/svs/igv/{}".format(self.outdir, sample_str)
        svict_igvinput.jobname = "generate-igvnav-input-svict"

        self.add(svict_igvinput)

    def configure_manta(self, normal_capture, cancer_capture):
        """
        Configure manta, to identify structural variants in sample

        :param normal_capture: A unique normal sample library capture
        :param cancer_capture: A unique cancer sample library capture
        """
        cancer_bam = self.get_capture_bam(cancer_capture, umi=False)
        normal_bam = self.get_capture_bam(normal_capture, umi=False)
        target_name = self.get_capture_name(cancer_capture.capture_kit_id)

        cancer_capture_str = compose_lib_capture_str(cancer_capture)
        normal_capture_str = compose_lib_capture_str(normal_capture)

        manta_sv = MantaSomaticSV()
        manta_sv.input_tumor = cancer_bam
        manta_sv.input_normal = normal_bam
        manta_sv.tumorid = cancer_capture_str
        manta_sv.normalid = normal_capture_str
        manta_sv.reference_sequence = self.refdata["reference_genome"]
        manta_sv.target_bed = self.refdata['targets'][target_name]['targets-bed-slopped20-gz']
        manta_sv.output_dir = "{}/svs/{}-{}-manta-somatic".format(self.outdir, normal_capture_str, cancer_capture_str)

        self.add(manta_sv)

    def configure_sv_calling(self, normal_capture, cancer_capture):
        """
        Configure Structural Variant Calling, to identify structural variants in sample

        :param normal_capture: A unique normal sample library capture
        :param cancer_capture: A unique cancer sample library capture
        """
        cancer_bam = self.get_capture_bam(cancer_capture, umi=False)
        normal_bam = self.get_capture_bam(normal_capture, umi=False)
        target_name = self.get_capture_name(cancer_capture.capture_kit_id)

        cancer_capture_str = compose_lib_capture_str(cancer_capture)
        normal_capture_str = compose_lib_capture_str(normal_capture)

        svaba = Svaba()
        svaba.scratch = self.scratch
        svaba.input_normal = normal_bam
        svaba.input_tumor = cancer_bam
        svaba.reference_sequence = self.refdata["bwaIndex"]
        svaba.threads = self.maxcores
        svaba.target_bed = self.refdata['targets'][target_name]['targets-bed-slopped20-gz']
        svaba.output_sample = "{}/svs/svaba/{}-{}".format(self.outdir, normal_capture_str, cancer_capture_str)
        svaba.output_somatic = "{}/svs/svaba/{}-{}.svaba.somatic.sv.vcf".format(self.outdir, normal_capture_str, cancer_capture_str)
        svaba.output_germline = "{}/svs/svaba/{}-{}.svaba.germline.sv.vcf".format(self.outdir, normal_capture_str, cancer_capture_str)

        self.add(svaba)

        annotate_svaba = {}
        svaba_igvinput = {}

        for vcftype in ['germline', 'somatic']:
            annotate_svaba[vcftype] = AnnotateSvaba()
            annotate_svaba[vcftype].input_vcf = svaba.output_somatic if vcftype == 'somatic' else svaba.output_germline
            annotate_svaba[vcftype].output_vcf = "{}/svs/svaba/{}-{}.svaba.{}.annotated.sv.vcf".format(self.outdir, 
                                                                                            normal_capture_str, 
                                                                                            cancer_capture_str,
                                                                                            vcftype)
            annotate_svaba[vcftype].jobname = "annotate-svtype-svaba-" + vcftype

            self.add(annotate_svaba[vcftype])

            svaba_igvinput[vcftype] = GenerateIGVNavInputSV()
            svaba_igvinput[vcftype].input_vcf = annotate_svaba[vcftype].output_vcf
            svaba_igvinput[vcftype].sdid = cancer_capture.sdid
            svaba_igvinput[vcftype].tool = 'svaba'
            svaba_igvinput[vcftype].vcftype = vcftype
            svaba_igvinput[vcftype].output = "{}/svs/igv/{}-{}".format(self.outdir, normal_capture_str, cancer_capture_str)
            svaba_igvinput[vcftype].output_mut = "{}/svs/igv/{}-{}_{}_svaba.mut".format(self.outdir, normal_capture_str, cancer_capture_str, vcftype)
            svaba_igvinput[vcftype].jobname = "generate-igvnav-input-svaba-" + vcftype

            self.add(svaba_igvinput[vcftype])


        lumpy = Lumpy()
        lumpy.scratch = self.scratch
        lumpy.input_normal = normal_bam
        lumpy.input_tumor = cancer_bam
        lumpy.normal_discordants = "{}/svs/lumpy/{}-discordants.bam".format(self.outdir, normal_capture_str)
        lumpy.tumor_discordants = "{}/svs/lumpy/{}-discordants.bam".format(self.outdir,  cancer_capture_str)
        lumpy.normal_splitters = "{}/svs/lumpy/{}-splitters.bam".format(self.outdir, normal_capture_str) 
        lumpy.tumor_splitters = "{}/svs/lumpy/{}-splitters.bam".format(self.outdir,  cancer_capture_str) 
        lumpy.output = "{}/svs/lumpy/{}-{}-lumpy.vcf".format(self.outdir, normal_capture_str, cancer_capture_str)
        lumpy.threads = self.maxcores

        self.add(lumpy)

        lumpy_igvinput = GenerateIGVNavInputSV()
        lumpy_igvinput.input_vcf = lumpy.output
        lumpy_igvinput.sdid = cancer_capture.sdid
        lumpy_igvinput.tool = 'lumpy'
        lumpy_igvinput.vcftype = 'somatic'
        lumpy_igvinput.output = "{}/svs/igv/{}-{}".format(self.outdir, normal_capture_str, cancer_capture_str)
        lumpy_igvinput.output_mut = "{}/svs/igv/{}-{}_lumpy_len500_SU24.mut".format(self.outdir, normal_capture_str, cancer_capture_str)
        lumpy_igvinput.jobname = "generate-igvnav-input-lumpy"

        self.add(lumpy_igvinput)

        annotate_sv = GenerateIGVNavInputSV()
        annotate_sv.input_vcf = "{}/svs/igv/".format(self.outdir)
        annotate_sv.input_svaba_germline = svaba_igvinput['germline'].output_mut
        annotate_sv.input_svaba_somatic = svaba_igvinput['somatic'].output_mut
        annotate_sv.input_lumpy = lumpy_igvinput.output_mut
        annotate_sv.input_svict_normal = "{}/svs/igv/{}_svict_SR8.mut".format(self.outdir, normal_capture_str)
        annotate_sv.input_svict_cancer = "{}/svs/igv/{}_svict_SR8.mut".format(self.outdir, cancer_capture_str)
        annotate_sv.input_svcaller_normal = "{}/svs/igv/{}_svcaller.mut".format(self.outdir, normal_capture_str)
        annotate_sv.input_svcaller_cancer = "{}/svs/igv/{}_svcaller.mut".format(self.outdir, cancer_capture_str)
        annotate_sv.input_gridss_normal = "{}/svs/igv/{}_pass_gridss.mut".format(self.outdir, normal_capture_str)
        annotate_sv.input_gridss_cancer = "{}/svs/igv/{}_pass_gridss.mut".format(self.outdir, cancer_capture_str)
        annotate_sv.annot_bed = self.refdata['genes_bed']
        annotate_sv.target_bed = self.refdata['targets'][target_name]['targets-bed-slopped20']
        annotate_sv.output = "{}/svs/igv/{}-{}-sv-annotated.txt".format(self.outdir, normal_capture_str, cancer_capture_str)
        annotate_sv.jobname = "annotate-strutural-variant"

        self.add(annotate_sv)


    def configure_liqbio_cna(self, normal_capture, cancer_capture):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)]
        tumor_results = self.capture_to_results[cancer_capture]

        # Configure the liqbio frankenplots:
        # NOTE: Get PureCN outputs:
        pureCN_outputs = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].pureCN_outputs

        normal_str = compose_lib_capture_str(normal_capture)
        cancer_str = compose_lib_capture_str(cancer_capture)

        liqbio_cna = LiqbioCNAPlot()
        liqbio_cna.input_tumor_cnr = self.capture_to_results[cancer_capture].cnr
        liqbio_cna.input_tumor_cns = self.capture_to_results[cancer_capture].cns
        liqbio_cna.input_normal_cnr = self.capture_to_results[normal_capture].cnr
        liqbio_cna.input_normal_cns = self.capture_to_results[normal_capture].cns
        liqbio_cna.input_het_snps_vcf = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vcf_addsample_output
        liqbio_cna.input_purecn_csv = pureCN_outputs["csv"]
        liqbio_cna.input_purecn_genes_csv = pureCN_outputs["genes_csv"]
        liqbio_cna.input_purecn_loh_csv = pureCN_outputs["loh_csv"]
        liqbio_cna.input_purecn_variants_csv = pureCN_outputs["variants_csv"]
        liqbio_cna.input_svcaller_T_DEL = self.capture_to_results[cancer_capture].svs["DEL"][1]
        liqbio_cna.input_svcaller_T_DUP = self.capture_to_results[cancer_capture].svs["DUP"][1]
        liqbio_cna.input_svcaller_T_INV = self.capture_to_results[cancer_capture].svs["INV"][1]
        liqbio_cna.input_svcaller_T_TRA = self.capture_to_results[cancer_capture].svs["TRA"][1]
        liqbio_cna.input_svcaller_N_DEL = self.capture_to_results[normal_capture].svs["DEL"][1]
        liqbio_cna.input_svcaller_N_DUP = self.capture_to_results[normal_capture].svs["DUP"][1]
        liqbio_cna.input_svcaller_N_INV = self.capture_to_results[normal_capture].svs["INV"][1]
        liqbio_cna.input_svcaller_N_TRA = self.capture_to_results[normal_capture].svs["TRA"][1]
        liqbio_cna.input_germline_mut_vcf = self.get_vepped_germline_vcf(normal_capture) # Vepped germline variants
        liqbio_cna.input_somatic_mut_vcf = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vepped_vcf
        liqbio_cna.output_plot_png = "{}/qc/{}-{}-liqbio-cna.png".format(self.outdir, normal_str, cancer_str)
        liqbio_cna.output_plot_png_normal = "{}/qc/{}-liqbio-cna.png".format(self.outdir, normal_str)
        liqbio_cna.output_cna_json = "{}/variants/{}-{}-liqbio-cna.json".format(self.outdir, normal_str, cancer_str)
        liqbio_cna.output_purity_json = "{}/qc/{}-{}-liqbio-purity.json".format(self.outdir, normal_str, cancer_str)
        self.add(liqbio_cna)

    def configure_panel_analysis_cancer_vs_normal_liqbio(self, normal_capture, cancer_capture, umi):
        capture_name = self.get_capture_name(cancer_capture.capture_kit_id)

        self.configure_sv_calling(normal_capture, cancer_capture)

        # self.configure_manta(normal_capture, cancer_capture)

        # PureCN doesn't require a aspecific targets/gene file any more, so it can be run always
        self.configure_purecn(normal_capture, cancer_capture, umi)
        self.configure_liqbio_cna(normal_capture, cancer_capture)

    def configure_umi_processing(self):
        # configure for UMI SNV calling pipeline
        #
        capture_to_barcodes = self.get_unique_capture_to_clinseq_barcodes()
        for unique_capture in capture_to_barcodes.keys():
            capture_kit = unique_capture.capture_kit_id
            for clinseq_barcode in capture_to_barcodes[unique_capture]:
                trimmed_fqfiles = fq_trimming(self,
                                  fq1_files=find_fastqs(clinseq_barcode, self.libdir)[0],
                                  fq2_files=find_fastqs(clinseq_barcode, self.libdir)[1],
                                  clinseq_barcode=clinseq_barcode,
                                  ref=self.refdata['bwaIndex'],
                                  outdir= "{}/bams/{}".format(self.outdir, capture_kit),
                                  maxcores=self.maxcores)
            
            bam_file = self.configure_fastq_to_bam(fq_files=trimmed_fqfiles, 
                                                    clinseq_barcode=clinseq_barcode, 
                                                    capture_kit=capture_kit)
            realigned_bam = self.configure_alignment_with_umi(bamfile=bam_file, 
                                                    clinseq_barcode=clinseq_barcode, 
                                                    capture_kit=capture_kit, jobname='1')
            consensus_reads = self.configure_consensus_reads_calling(bam=realigned_bam, 
                                                    clinseq_barcode=clinseq_barcode,
                                                    capture_kit=capture_kit)
            realigned_bam2 = self.configure_alignment_with_umi(bamfile=consensus_reads, 
                                                    clinseq_barcode=clinseq_barcode, 
                                                    capture_kit=capture_kit, jobname='2')
            filtered_bam = self.configure_consensus_read_filter(bam=realigned_bam2 ,
                                                    clinseq_barcode=clinseq_barcode,
                                                    capture_kit=capture_kit)
            clip_overlap_bam = self.configure_clip_overlapping(bam=filtered_bam,
                                                    clinseq_barcode=clinseq_barcode,
                                                    capture_kit=capture_kit)
            mark_dups_bam = self.configure_markdups(bamfile=realigned_bam, unique_capture=unique_capture)

            self.set_capture_bam(unique_capture, clip_overlap_bam, self.umi)

    def configure_alignment_with_umi(self, bamfile, clinseq_barcode, capture_kit, jobname):
        # Map the reads with bwa and merge with the UMI tags (picard SamToFastq | bwa mem | picard MergeBamAlignment)
        align_unmap_bam = AlignUnmappedBam()
        align_unmap_bam.input_bam = bamfile
        align_unmap_bam.reference_genome = self.refdata['bwaIndex']
        align_unmap_bam.threads = self.maxcores
        align_unmap_bam.scratch = self.scratch
        align_unmap_bam.is_intermediate = True
        align_unmap_bam.output_bam = "{}/bams/{}/{}.mapped-{}.bam".format(self.outdir, capture_kit, clinseq_barcode, jobname)
        align_unmap_bam.jobname = "alignment-of-unmapped-bam-"+ jobname + '-' + clinseq_barcode
        self.add(align_unmap_bam)

        targets = self.get_capture_name(capture_kit)

        realignment = Realignment()
        realignment.input_bam = align_unmap_bam.output_bam
        realignment.scratch = self.scratch
        realignment.output_bam = "{}/bams/{}/{}.realigned-{}.bam".format(self.outdir, capture_kit, clinseq_barcode, jobname)
        realignment.reference_genome = self.refdata['reference_genome']
        realignment.target_region = self.refdata['targets'][targets]['targets-bed-slopped20']
        realignment.known_indel1 = self.refdata['1KG']
        realignment.known_indel2 = self.refdata['Mills_and_1KG_gold_standard']
        realignment.target_intervals = "{}/bams/{}/{}.intervals".format(self.outdir, capture_kit, clinseq_barcode)
        realignment.is_intermediate = True
        realignment.jobname = "realignment-" + jobname + '-' + clinseq_barcode
        self.add(realignment)

        return realignment.output_bam

    def configure_fastq_to_bam(self, fq_files, clinseq_barcode, capture_kit):
        # Extract UMIs from trimmed fastq and store in RX tag of unmapped bam (fgbio FastqToBam)
        
        library = parse_prep_id(clinseq_barcode)
        sample = compose_sample_str(extract_unique_capture(clinseq_barcode))

        fastq_to_bam = FastqToBam()
        fastq_to_bam.input_fastq1 = fq_files[0]
        fastq_to_bam.input_fastq2 = fq_files[1]
        fastq_to_bam.sample = sample
        fastq_to_bam.library = library
        fastq_to_bam.scratch = self.scratch
        fastq_to_bam.is_intermediate = True
        fastq_to_bam.output_bam = "{}/bams/{}/{}.unmapped.bam".format(self.outdir, capture_kit, clinseq_barcode)
        fastq_to_bam.jobname = "fastq-to-bam" + '-' + clinseq_barcode
        self.add(fastq_to_bam)

        return fastq_to_bam.output_bam

    def configure_consensus_reads_calling(self, bam,  clinseq_barcode, capture_kit):

        group_reads = GroupReadsByUmi()
        group_reads.input_bam = bam
        group_reads.scratch = self.scratch
        group_reads.histogram = "{}/bams/{}/{}.grouped.bam.fs.txt".format(self.outdir, capture_kit, clinseq_barcode)
        group_reads.output_bam = "{}/bams/{}/{}.grouped.bam".format(self.outdir, capture_kit, clinseq_barcode)
        group_reads.is_intermediate = True
        group_reads.jobname = "group-reads-by-umi" + '-' + clinseq_barcode
        self.add(group_reads)

        call_consensus_reads = CallDuplexConsensusReads()
        call_consensus_reads.input_bam = group_reads.output_bam
        call_consensus_reads.scratch = self.scratch
        call_consensus_reads.threads = self.maxcores
        call_consensus_reads.output_bam = "{}/bams/{}/{}.consensus.bam".format(self.outdir, capture_kit, clinseq_barcode)
        call_consensus_reads.is_intermediate = True
        call_consensus_reads.jobname = "call-duplex-consensus-reads" + '-' + clinseq_barcode
        self.add(call_consensus_reads)

        return call_consensus_reads.output_bam

    def configure_consensus_read_filter(self, bam, clinseq_barcode, capture_kit):

        filter_con_reads = FilterConsensusReads()
        filter_con_reads.input_bam = bam
        filter_con_reads.scratch = self.scratch
        filter_con_reads.reference_genome = self.refdata['reference_genome']
        filter_con_reads.output_bam = "{}/bams/{}/{}.consensus.filtered.bam".format(self.outdir, capture_kit, clinseq_barcode)
        filter_con_reads.jobname = "filter-consensus-reads-{}".format(clinseq_barcode)
        self.add(filter_con_reads)

        return filter_con_reads.output_bam

    def configure_clip_overlapping(self, bam, clinseq_barcode, capture_kit):

        clip_overlap_reads = ClipBam()
        clip_overlap_reads.input_bam = bam
        clip_overlap_reads.scratch = self.scratch
        clip_overlap_reads.reference_genome = self.refdata['reference_genome']
        clip_overlap_reads.output_bam = "{}/bams/{}/{}.clip.overlapped.bam".format(self.outdir, capture_kit, clinseq_barcode)
        clip_overlap_reads.metrics_txt = "{}/qc/{}-clip_overlap_metrix.txt".format(self.outdir, clinseq_barcode)
        clip_overlap_reads.jobname = "clip-overlap-reads-{}".format(clinseq_barcode)
        self.add(clip_overlap_reads)

        return clip_overlap_reads.output_bam




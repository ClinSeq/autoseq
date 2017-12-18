from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.cnvcalling import LiqbioCNAPlot
from autoseq.util.clinseq_barcode import *
from autoseq.tools.structuralvariants import Svcaller, Sveffect

__author__ = 'thowhi'


class LiqBioPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, job_params, outdir, libdir, maxcores=1, scratch="/scratch/tmp/tmp",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, job_params, outdir, libdir,
                                 maxcores, scratch, **kwargs)

        # Set the min alt frac value:
        self.default_job_params["vardict-min-alt-frac"] = 0.01
        self.default_job_params["vardict-min-num-reads"] = None
        self.default_job_params["vep-additional-options"] = " --pick --filter_common "

        # Remove clinseq barcodes for which data is not available:
        self.check_sampledata()

        # Configure alignment and merging of fastq data for all clinseq barcodes:
        self.configure_align_and_merge()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure liqbio-specific panel analyses:
        self.configure_panel_analyses_liqbio()

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

    def configure_single_capture_analysis_liqbio(self, unique_capture):
        input_bam = self.get_capture_bam(unique_capture)
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
        sveffect.input_del_gtf = self.capture_to_results[unique_capture]["DEL"][1]
        sveffect.input_dup_gtf = self.capture_to_results[unique_capture]["DUP"][1]
        sveffect.input_inv_gtf = self.capture_to_results[unique_capture]["INV"][1]
        sveffect.input_tra_gtf = self.capture_to_results[unique_capture]["TRA"][1]
        sveffect.ts_regions = self.refdata["ts_regions"]
        sveffect.ar_regions = self.refdata["ar_regions"]
        sveffect.fusion_regions = self.refdata["fusion_regions"]
        sveffect.output_combined_bed = "{}/svs/{}_combined.bed".format(self.outdir, sample_str)
        sveffect.output_effects_json = "{}/svs/{}_effects.json".format(self.outdir, sample_str)

        self.add(sveffect)

        self.set_capture_sveffect(unique_capture, sveffect.output_effects_json)

    def configure_panel_analyses_liqbio(self):
        # Configure liqbio analyses to be run on all unique panel captures individually:
        for unique_capture in self.get_mapped_captures_no_wgs():
            self.configure_single_capture_analysis_liqbio(unique_capture)

        # Configure a liqbio analyses for each normal-cancer pairing:
        for normal_capture in self.get_mapped_captures_normal():
            for cancer_capture in self.get_mapped_captures_cancer():
                self.configure_panel_analysis_cancer_vs_normal_liqbio(
                    normal_capture, cancer_capture)

    def configure_liqbio_cna(self, normal_capture, cancer_capture):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)]
        tumor_results = self.capture_to_results[cancer_capture]

        # Configure the liqbio frankenplots:
        # NOTE: Get PureCN outputs:
        pureCN_outputs = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].pureCN_outputs

        cancer_str = compose_lib_capture_str(cancer_capture)

        liqbio_cna = LiqbioCNAPlot()
        liqbio_cna.tumor_cnr = self.capture_to_results[cancer_capture].cnr
        liqbio_cna.tumor_cns = self.capture_to_results[cancer_capture].cns
        liqbio_cna.normal_cnr = self.capture_to_results[normal_capture].cnr
        liqbio_cna.normal_cns = self.capture_to_results[normal_capture].cns
        liqbio_cna.het_snps_vcf = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vcf_addsample_output
        liqbio_cna.purecn_csv = pureCN_outputs["csv"]
        liqbio_cna.purecn_genes_csv = pureCN_outputs["genes_csv"]
        liqbio_cna.purecn_loh_csv = pureCN_outputs["loh_csv"]
        liqbio_cna.purecn_variants_csv = pureCN_outputs["variants_csv"]
        liqbio_cna.svcaller_T_DEL = self.capture_to_results[cancer_capture]["DEL"][1]
        liqbio_cna.svcaller_T_DUP = self.capture_to_results[cancer_capture]["DUP"][1]
        liqbio_cna.svcaller_T_INV = self.capture_to_results[cancer_capture]["INV"][1]
        liqbio_cna.svcaller_T_TRA = self.capture_to_results[cancer_capture]["TRA"][1]
        liqbio_cna.svcaller_N_DEL = self.capture_to_results[normal_capture]["DEL"][1]
        liqbio_cna.svcaller_N_DUP = self.capture_to_results[normal_capture]["DUP"][1]
        liqbio_cna.svcaller_N_INV = self.capture_to_results[normal_capture]["INV"][1]
        liqbio_cna.svcaller_N_TRA = self.capture_to_results[normal_capture]["TRA"][1]
        liqbio_cna.germline_mut_vcf = self.normal_capture_to_vcf[normal_capture]
        liqbio_cna.somatic_mut_vcf = self.normal_cancer_pair_to_results[(normal_capture, cancer_capture)].vepped_vcf
        liqbio_cna.plot_png = "{}/cnv/{}-{}_combined.bed".format(self.outdir, normal_capture, cancer_capture)
        liqbio_cna.output_cna = # XXX FIGURE THIS OUT: LOOK AT ALASCCA OUTPUTS.
        liqbio_cna.output_purity = # XXX CONTINUE HERE: FILL THESE IN

        return liqbio_cna.output_cna, liqbio_cna.output_purity

    def configure_panel_analysis_cancer_vs_normal_liqbio(self, normal_capture, cancer_capture):
        cancer_capture_str = compose_lib_capture_str(cancer_capture)

        # XXX CONFIGURE SVCALLER, AND STORE THE OUTPUT FILES SOMEHOW -> 

        if self.refdata['targets'][cancer_capture_str]['purecn_targets']:
            self.configure_purecn(normal_capture, cancer_capture)

            self.configure_liqbio_cna(normal_capture, cancer_capture)
import logging

from autoseq.pipeline.clinseq import ClinseqPipeline
from autoseq.tools.cnvcalling import AlasccaCNAPlot
from autoseq.tools.reports import CompileMetadata, CompileAlasccaGenomicJson, WriteAlasccaReport
from autoseq.pipeline.clinseq import compose_lib_capture_str
from autoseq.util.report_type import only_alascca_class_report

__author__ = 'dankle'


class AlasccaPipeline(ClinseqPipeline):
    def __init__(self, sampledata, refdata, job_params, outdir, libdir, maxcores=1, scratch="/scratch/tmp/tmp/",
                 referral_db_conf="tests/referrals/referral-db-config.json",
                 addresses="tests/referrals/addresses.csv",
                 **kwargs):
        ClinseqPipeline.__init__(self, sampledata, refdata, job_params, outdir, libdir,
                                 maxcores, scratch, **kwargs)

        self.referral_db_conf = referral_db_conf
        self.addresses = addresses
        self.default_job_params["vardict-min-num-reads"] = 6
        self.default_job_params["create_alascca_report"] = True

        # Check to ensure that the sample data is valid for an ALASCCA analysis:
        self.validate_sample_data_for_alascca()

        # Remove sample capture items for which data is not available:
        self.check_sampledata()

        # Configure alignment and merging of fastq data for all clinseq barcodes:
        self.configure_align_and_merge()

        # Configure all panel analyses:
        self.configure_panel_analyses()

        # Configure QC of all panel data:
        self.configure_all_panel_qcs()

        # Configure ALASCCA report generation:
        self.configure_alascca_specific_analysis()

        # Configure fastq QCs:
        self.configure_fastq_qcs()

        # Configure MultiQC:
        self.configure_multi_qc()

    def validate_sample_data_for_alascca(self):
        """
        Checks validity of the sample data. Raises a ValueError if the sampledata dictionary does
        not fit into the expected ALASCCA analysis pipeline limitations.
        """

        # FIXME: Implement this
        pass

    def get_normal_and_tumor_captures(self):
        """
        Retrieves the unique normal and tumor capture identifiers for this ALASCCA analysis.
        :return: (normal_capture, tumor_capture) tuple, denoting those unique library captures.
        """

        # There must be exactly one tumor and exactly one normal for this to be valid:
        if len(self.get_mapped_captures_cancer()) != 1 or \
           len(self.get_mapped_captures_normal()) != 1:
            raise ValueError("Invalid pipeline state for configuration of ALASCCA CNA.")

        normal_capture = self.get_mapped_captures_normal()[0]
        tumor_capture = self.get_mapped_captures_cancer()[0]

        return normal_capture, tumor_capture

    def configure_alascca_cna(self, normal_capture, tumor_capture):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, tumor_capture)]
        tumor_results = self.capture_to_results[tumor_capture]

        alascca_cna = AlasccaCNAPlot()
        alascca_cna.input_somatic_vcf = tumor_vs_normal_results.vepped_vcf
        alascca_cna.input_germline_vcf = tumor_vs_normal_results.vcf_addsample_output
        alascca_cna.input_cnr = tumor_results.cnr
        alascca_cna.input_cns = tumor_results.cns
        alascca_cna.chrsizes = self.refdata['chrsizes']

        tumor_str = compose_lib_capture_str(tumor_capture)

        alascca_cna.output_cna = "{}/variants/{}-alascca-cna.json".format(
            self.outdir, tumor_str)
        alascca_cna.output_purity = "{}/variants/{}-alascca-purity.json".format(
            self.outdir, tumor_str)
        alascca_cna.output_png = "{}/qc/{}-alascca-cna.png".format(
            self.outdir, tumor_str)
        alascca_cna.jobname = "alascca-cna/{}".format(tumor_str)
        self.add(alascca_cna)

        return alascca_cna.output_cna, alascca_cna.output_purity

    def configure_compile_metadata(self, normal_capture, tumor_capture):
        blood_barcode = normal_capture.sample_id
        tumor_barcode = tumor_capture.sample_id
        metadata_json = "{}/report/{}-{}.metadata.json".format(self.outdir, blood_barcode, tumor_barcode)
        compile_metadata_json = CompileMetadata(
            self.referral_db_conf, blood_barcode, tumor_barcode, metadata_json, self.addresses)
        compile_metadata_json.jobname = "compile-metadata/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_metadata_json)

        return compile_metadata_json.output_json

    def configure_compile_genomic_json(self, normal_capture, tumor_capture,
                                       alascca_cna_output, alascca_cna_purity_call):
        tumor_vs_normal_results = self.normal_cancer_pair_to_results[(normal_capture, tumor_capture)]
        tumor_results = self.capture_to_results[tumor_capture]
        normal_results = self.capture_to_results[normal_capture]

        blood_barcode = normal_capture.sample_id
        tumor_barcode = tumor_capture.sample_id

        genomic_json = "{}/report/{}-{}.genomic.json".format(
            self.outdir, blood_barcode, tumor_barcode)

        compile_genomic_json = CompileAlasccaGenomicJson(
            input_somatic_vcf=tumor_vs_normal_results.vepped_vcf,
            input_cn_calls=alascca_cna_output,
            input_msisensor=tumor_vs_normal_results.msi_output,
            input_purity_qc=alascca_cna_purity_call,
            input_contam_qc=tumor_vs_normal_results.cancer_contam_call,
            input_tcov_qc=tumor_results.cov_qc_call,
            input_ncov_qc=normal_results.cov_qc_call,
            output_json=genomic_json)

        compile_genomic_json.jobname = "compile-genomic/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(compile_genomic_json)

        return compile_genomic_json.output_json

    def configure_write_alascca_report(
             self, normal_capture, tumor_capture, metadata_json, genomic_json):
        blood_barcode = normal_capture.sample_id
        tumor_barcode = tumor_capture.sample_id
        only_alascca = only_alascca_class_report(blood_barcode, tumor_barcode, self.referral_db_conf)

        pdf = "{}/report/AlasccaReport-{}-{}.pdf".format(self.outdir, blood_barcode, tumor_barcode)
        write_alascca_pdf = WriteAlasccaReport(
            genomic_json, metadata_json, pdf, only_alascca)
        write_alascca_pdf.jobname = "writeAlasccaPdf/{}-{}".format(tumor_barcode, blood_barcode)
        self.add(write_alascca_pdf)

    def configure_alascca_report_generation(self, normal_capture, tumor_capture,
                                            alascca_cna_output, alascca_cna_purity_call):
        """
        Configure the generation of the ALASCCA report for this pipeline instance.
        """

        # Always create genomic json
        genomic_json = self.configure_compile_genomic_json(
            normal_capture, tumor_capture, alascca_cna_output, alascca_cna_purity_call)

        # Create metadata json and final pdf only if supposed to
        create_alascca_report = self.get_job_param("create_alascca_report")
        if create_alascca_report:
            metadata_json = self.configure_compile_metadata(normal_capture, tumor_capture)
            self.configure_write_alascca_report(normal_capture, tumor_capture,
                                                metadata_json, genomic_json)

    def configure_alascca_specific_analysis(self):
        """
        Configure the Jobs for specific to the ALASCCA pipeline. The panel analyses
        and panel QC analyses must be configured before this can be configured.
        """

        normal_capture, tumor_capture = self.get_normal_and_tumor_captures()
        alascca_cna_ouput, alascca_cna_purity = \
            self.configure_alascca_cna(normal_capture, tumor_capture)
        self.configure_alascca_report_generation(
            normal_capture, tumor_capture, alascca_cna_ouput, alascca_cna_purity)

    def get_coverage_bed(self, targets):
        """
        Retrieve the targets bed file to use for calculating coverage, given the specified
        targets name. Uses the alascca bed file if available.

        :param targets: Target capture name
        :return: bed file name
        """
        if 'alascca_targets' in self.refdata['targets']:
            return self.refdata['targets']['alascca_targets']['targets-bed-slopped20']
        else:
            return self.refdata['targets'][targets]['targets-bed-slopped20']

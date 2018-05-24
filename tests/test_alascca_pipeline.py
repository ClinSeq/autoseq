import unittest
from mock import patch
from autoseq.pipeline.alascca import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestAlascca(unittest.TestCase):
    def setUp(self):
        self.sample_data_invalid = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
        }

        self.sample_data_valid = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1"],
            "CFDNA": []
        }

        self.ref_data = {
            "bwaIndex": "bwa/test-genome-masked.fasta",
            "chrsizes": "genome/test-genome-masked.chrsizes.txt",
            "clinvar": "variants/clinvar_20160203.vcf.gz",
            "cosmic": "variants/CosmicCodingMuts_v71.vcf.gz",
            "dbSNP": "variants/dbsnp142-germline-only.vcf.gz",
            "exac": "variants/ExAC.r0.3.1.sites.vep.vcf.gz",
            "icgc": "variants/icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz",
            "reference_dict": "genome/test-genome-masked.dict",
            "reference_genome": "genome/test-genome-masked.fasta",
            "swegene_common": "variants/swegen_common.vcf.gz",
            "targets": {
                "test-regions": {
                    "cnvkit-ref": {
                        "THRUPLEX_PLASMASEQ": {
                            "CFDNA": "intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnn", 
                            "N": "intervals/targets/progression.THRUPLEX_PLASMASEQ.N.cnn"
                        }
                    },
                    "cnvkit-fix": {
                        "THRUPLEX_PLASMASEQ": {
                            "CFDNA": "intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnvkit-fix.tsv" 
                        }
                    },
                    "msisites": "intervals/targets/test-regions.msisites.tsv",
                    "targets-bed-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.bed",
                    "targets-interval_list": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "targets-interval_list-slopped20": "intervals/targets/test-regions-GRCh37.slopped20.interval_list",
                    "blacklist-bed": None,
                }
            },
            "contest_vcfs": {
                "test-regions": "test_contest.vcf"
            },
            "vep_dir": "dummy_vep_dir"
        }

        self.test_tumor_capture = UniqueCapture("AL", "P-NA12877", "T", "03098849", "TD", "TT")
        self.test_normal_capture = UniqueCapture("AL", "P-NA12877", "N", "03098121", "TD", "TT")

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    @patch('autoseq.pipeline.alascca.AlasccaPipeline.configure_multi_qc')
    def test_constructor_valid(self, mock_configure_multi_qc, mock_find_fastqs,
                               mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        self.assertTrue(mock_configure_multi_qc.called)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_constructor_invalid(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        self.assertRaises(ValueError, lambda: AlasccaPipeline(self.sample_data_invalid, self.ref_data,
                                                              {}, "/tmp", "/nfs/LIQBIO/INBOX/exomes"))

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_get_normal_and_tumor_captures(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        normal_capture, tumor_capture = test_alascca_pipeline.get_normal_and_tumor_captures()
        self.assertEquals(normal_capture, self.test_normal_capture)
        self.assertEquals(tumor_capture, self.test_tumor_capture)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_alascca_cna(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        output_cna, output_purity = test_alascca_pipeline.configure_alascca_cna(self.test_normal_capture,
                                                                                self.test_tumor_capture)
        self.assertEquals(output_cna.split(".")[-1], "json")
        self.assertEquals(output_purity.split(".")[-1], "json")

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_compile_metadata(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        metadata_json = test_alascca_pipeline.configure_compile_metadata(self.test_normal_capture,
                                                                         self.test_tumor_capture)
        self.assertEquals(metadata_json.split(".")[-1], "json")

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_compile_genomic_json(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        genomic_json = test_alascca_pipeline.configure_compile_genomic_json(self.test_normal_capture,
                                                                            self.test_tumor_capture,
                                                                            "dummy.txt", "dummy.json")
        self.assertEquals(genomic_json.split(".")[-1], "json")

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_write_alascca_report(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        num_jobs_before_call = len(test_alascca_pipeline.graph.nodes())
        genomic_json = test_alascca_pipeline.configure_write_alascca_report(self.test_normal_capture,
                                                                            self.test_tumor_capture,
                                                                            "dummy.txt", "dummy.json")
        num_jobs_after_call = len(test_alascca_pipeline.graph.nodes())
        self.assertEquals(num_jobs_before_call + 1, num_jobs_after_call)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_alascca_report_generation(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        num_jobs_before_call = len(test_alascca_pipeline.graph.nodes())
        test_alascca_pipeline.configure_alascca_report_generation(self.test_normal_capture,
                                                                  self.test_tumor_capture,
                                                                  "dummy.txt", "dummy.json")
        num_jobs_after_call = len(test_alascca_pipeline.graph.nodes())
        self.assertEquals(num_jobs_before_call + 3, num_jobs_after_call)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_alascca_report_generation_no_report(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {"create_alascca_report": False},
                                                "/tmp", "/nfs/LIQBIO/INBOX/exomes")
        num_jobs_before_call = len(test_alascca_pipeline.graph.nodes())
        test_alascca_pipeline.configure_alascca_report_generation(self.test_normal_capture,
                                                                  self.test_tumor_capture,
                                                                  "dummy.txt", "dummy.json")
        num_jobs_after_call = len(test_alascca_pipeline.graph.nodes())
        self.assertEquals(num_jobs_before_call + 1, num_jobs_after_call)

    @patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
    @patch('autoseq.pipeline.clinseq.find_fastqs')
    def test_configure_alascca_specific_analysis(self, mock_find_fastqs, mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_alascca_pipeline = AlasccaPipeline(self.sample_data_valid, self.ref_data, {}, "/tmp",
                                                "/nfs/LIQBIO/INBOX/exomes")
        num_jobs_before_call = len(test_alascca_pipeline.graph.nodes())
        test_alascca_pipeline.configure_alascca_specific_analysis()
        num_jobs_after_call = len(test_alascca_pipeline.graph.nodes())
        self.assertEquals(num_jobs_before_call + 4, num_jobs_after_call)

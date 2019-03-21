import unittest
from mock import patch
from autoseq.pipeline.liqbio import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestLiqbio(unittest.TestCase):
    def setUp(self):
        self.sample_data = {
            "sdid": "P-NA12877",
            "T": ["AL-P-NA12877-T-03098849-TD1-TT1", "AL-P-NA12877-T-03098849-TD1-WGS"],
            "N": ["AL-P-NA12877-N-03098121-TD1-TT1", "AL-P-NA12877-N-03098121-TD1-WGS"],
            "CFDNA": ["LB-P-NA12877-CFDNA-03098850-TD1-TT1", "LB-P-NA12877-CFDNA-03098850-TD1-TT2",
                      "LB-P-NA12877-CFDNA-03098850-TD1-WGS"]
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
            "ar_regions": "intervals/ar_regions.bed",
            "ts_regions": "intervals/ts_regions.bed",
            "fusion_regions": "intervals/fusion_regions.bed",
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
                    "purecn_targets": "intervals/targets/purecn.bed",
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
    @patch('autoseq.pipeline.liqbio.LiqBioPipeline.configure_multi_qc')
    def test_constructor_valid(self, mock_configure_multi_qc, mock_find_fastqs,
                               mock_data_available_for_clinseq_barcode):
        mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
        mock_data_available_for_clinseq_barcode.return_value = True
        test_liqbio_pipeline = LiqBioPipeline(self.sample_data, self.ref_data, {}, "/tmp",
                                              "/nfs/LIQBIO/INBOX/exomes", 'FALSE')
        self.assertTrue(mock_configure_multi_qc.called)

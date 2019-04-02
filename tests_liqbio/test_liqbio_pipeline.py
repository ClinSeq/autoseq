import unittest
from mock import patch
from autoseq.pipeline.liqbio import *
from autoseq.util.clinseq_barcode import UniqueCapture


class TestLiqbio(unittest.TestCase):
	def setUp(self):
		self.sample_data = {
			"sdid": "P-00202345",
			"N": [
				"LB-P-00202345-N-03277090-TP20190201-CP20190204"
			], 
			"T": [], 
			"CFDNA": [
				"LB-P-00202345-CFDNA-03277089-TP20190201-CP20190204", 
				"LB-P-00202345-CFDNA-03277089-TP20190201-CM20190204"
			]
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
			"1KG":"/nfs/ALASCCA/autoseq-genome/variants/1000G_phase1.indels.b37.vcf.gz", 
			"Mills_and_1KG_gold_standard": "/nfs/ALASCCA/autoseq-genome/variants/Mills_and_1000G_gold_standard.indels.b37.vcf.gz" ,
			"brca_exchange": "/nfs/ALASCCA/autoseq-genome/variants/BrcaExchangeClinvar_15Jan2019_v26_hg19.vcf.gz",
			"oncokb": "/nfs/ALASCCA/autoseq-genome/variants/OncoKB_6Mar19_v1.9.txt",
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
				},
				"progression": {
					"blacklist-bed": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.blacklist.bed", 
					"cnvkit-fix": {
						"THRUPLEX_PLASMASEQ": {
							"CFDNA": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnvkit-fix.tsv"
						}
					}, 
					"cnvkit-ref": {
						"KAPA_HYPERPREP": {
							"N": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.KAPA_HYPERPREP.N.cnn", 
							"T": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.KAPA_HYPERPREP.T.cnn"
						}, 
						"THRUPLEX_PLASMASEQ": {
							"CFDNA": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.THRUPLEX_PLASMASEQ.CFDNA.cnn", 
							"N": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.THRUPLEX_PLASMASEQ.N.cnn"
						}
					}, 
					"msings-baseline": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.msings.baseline", 
					"msings-bed": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.msings.bed", 
					"msings-msi_intervals": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.msings.msi_intervals", 
					"msisites": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.slopped20.msisites.tsv", 
					"purecn_targets": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.purecn.txt", 
					"targets-bed-slopped20": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.slopped20.bed.gz", 
					"targets-interval_list": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.interval_list", 
					"targets-interval_list-slopped20": "/nfs/ALASCCA/autoseq-genome/intervals/targets/progression.slopped20.interval_list"
				},
				"monitor": {
					"blacklist-bed": "/nfs/ALASCCA/autoseq-genome/intervals/targets/monitor.blacklist.bed", 
					"msings-baseline": None, 
					"msings-bed": None, 
					"msings-msi_intervals": None, 
					"msisites": "/nfs/ALASCCA/autoseq-genome/intervals/targets/monitor.slopped20.msisites.tsv", 
					"purecn_targets": None, 
					"targets-bed-slopped20": "/nfs/ALASCCA/autoseq-genome/intervals/targets/monitor.slopped20.bed.gz", 
					"targets-interval_list": "/nfs/ALASCCA/autoseq-genome/intervals/targets/monitor.interval_list", 
					"targets-interval_list-slopped20": "/nfs/ALASCCA/autoseq-genome/intervals/targets/monitor.slopped20.interval_list"
				}
			},
			"contest_vcfs": {
				"test-regions": "test_contest.vcf"
			},
			"vep_dir": "dummy_vep_dir"
		}
					
		self.test_tumor_capture = UniqueCapture("LB", "P-00202345", "CFDNA", "03277089", "TP", "CP")
		self.test_normal_capture = UniqueCapture("LB", "P-00202345", "N", "03277090", "TP", "CP")
		self.test_liqbio_pipeline = LiqBioPipeline(self.sample_data, self.ref_data, {}, "/tmp",
												"/media/clinseq/disk4/PROBIO/", False)

	@patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
	@patch('autoseq.pipeline.clinseq.find_fastqs')
	@patch('autoseq.pipeline.liqbio.LiqBioPipeline.configure_multi_qc')
	def test_constructor_valid(self, mock_configure_multi_qc, mock_find_fastqs,
							   mock_data_available_for_clinseq_barcode):
		mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
		mock_data_available_for_clinseq_barcode.return_value = True
		test_liqbio_pipeline_without_umi = LiqBioPipeline(self.sample_data, self.ref_data, {}, "/tmp/liqbio-test/",
											  "/media/clinseq/disk4/PROBIO/", False)
		self.assertTrue(mock_configure_multi_qc.called)

	@patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
	@patch('autoseq.pipeline.clinseq.find_fastqs')
	@patch('autoseq.pipeline.liqbio.LiqBioPipeline.configure_multi_qc')
	def test_constructor_valid_with_umi(self, mock_configure_multi_qc, mock_find_fastqs,
							   mock_data_available_for_clinseq_barcode):
		mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
		mock_data_available_for_clinseq_barcode.return_value = True
		test_liqbio_pipeline_with_umi = LiqBioPipeline(self.sample_data, self.ref_data, {}, "/tmp/liqbio-test/",
											  "/media/clinseq/disk4/PROBIO/", True)
		self.assertTrue(mock_configure_multi_qc.called)

	def test_configure_single_capture_analysis_liqbio(self):
		num_of_nodes_before = len(self.test_liqbio_pipeline.graph.nodes())
		self.test_liqbio_pipeline.configure_single_capture_analysis_liqbio(self.test_tumor_capture)
		num_of_nodes_after = len(self.test_liqbio_pipeline.graph.nodes())
		self.assertEquals(num_of_nodes_before+5, num_of_nodes_after)

	def  test_configure_panel_analyses_liqbio(self):
		num_of_nodes_before = len(self.test_liqbio_pipeline.graph.nodes())
		self.test_liqbio_pipeline.configure_panel_analyses_liqbio()
		num_of_nodes_after = len(self.test_liqbio_pipeline.graph.nodes())
		self.assertEquals(self.test_liqbio_pipeline.capture_to_results[self.test_tumor_capture].sv_effects.split('.')[-1], 'json')
		self.assertEquals(num_of_nodes_before+17, num_of_nodes_after)

	@patch('autoseq.pipeline.clinseq.data_available_for_clinseq_barcode')
	@patch('autoseq.pipeline.clinseq.find_fastqs')
	@patch('autoseq.tools.alignment.fq_trimming')
	def test_configure_umi_processing(self, mock_fq_trimming, mock_find_fastqs,
							   mock_data_available_for_clinseq_barcode):
		mock_find_fastqs.return_value = ["test1.fq.gz", "test2.fq.gz"]
		mock_data_available_for_clinseq_barcode.return_value = True
		mock_fq_trimming.return_value = ["trim_test1.fq.gz", "trim_test2.fq.gz"]
		test_liqbio_pipeline_with_umi = LiqBioPipeline(self.sample_data, self.ref_data, {}, "/tmp/liqbio-test/",
											  "/media/clinseq/disk4/PROBIO/", True)
		test_liqbio_pipeline_with_umi.configure_umi_processing()
		self.assertEquals(test_liqbio_pipeline_with_umi.capture_to_results[self.test_tumor_capture].umi_bamfile.split('.')[-1], 'bam')
		self.assertEquals(test_liqbio_pipeline_with_umi.capture_to_results[self.test_tumor_capture].merged_bamfile.split('.')[-1], 'bam')

	
		
		

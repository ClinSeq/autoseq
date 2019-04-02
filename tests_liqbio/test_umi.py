import unittest
from autoseq.tools.umi import *

class TestUMIcmds(unittest.TestCase):
	"""docstring for TestUMunittest.TestCase"""
	def setup(self):
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

	def test_fastqtobam(self):
		fastq_to_bam = FastqToBam()
		fastq_to_bam.input_fastq1 = 'input_1.fq'
		fastq_to_bam.input_fastq2 = 'input_2.fq'
		fastq_to_bam.output_bam = 'output.bam'
		fastq_to_bam.sample = 'P-SAMPLE'
		fastq_to_bam.library = 'TT'
		cmd = fastq_to_bam.command()
		self.assertIn('input_1.fq', cmd)
		self.assertIn('input_2.fq', cmd)
		self.assertIn('output.bam', cmd)
		self.assertIn('P-SAMPLE', cmd)
		self.assertIn('TT', cmd)

	def test_align_unmap_bam(self):
		align_unmap_bam = AlignUnmappedBam()
		align_unmap_bam.input_bam = 'input.bam'
		align_unmap_bam.reference_genome = 'human_g1k_v37_decoy.fasta'
		align_unmap_bam.output_bam = 'output.bam'
		cmd = align_unmap_bam.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('human_g1k_v37_decoy.fasta', cmd)
		self.assertIn('output.bam', cmd)

	def test_group_readsby_umi(self):
		group_reads = GroupReadsByUmi()
		group_reads.input_bam = 'input.bam'
		group_reads.output_bam = 'output.bam'
		group_reads.output_histogram = 'grouped.fs.txt'
		cmd = group_reads.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('output.bam', cmd)
		self.assertIn('grouped.fs.txt', cmd)

	def test_call_consensus_reads(self):
		call_consensus_reads = CallDuplexConsensusReads()
		call_consensus_reads.input_bam = 'input.bam'
		call_consensus_reads.output_bam = 'output.bam'
		cmd = call_consensus_reads.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('output.bam', cmd)

	def test_filter_consensus_reads(self):
		filter_consensus_reads = FilterConsensusReads()
		filter_consensus_reads.input_bam = 'input.bam'
		filter_consensus_reads.output_bam = 'output.bam'
		filter_consensus_reads.reference_genome = 'human_g1k_v37_decoy.fasta'
		cmd = filter_consensus_reads.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('output.bam', cmd)
		self.assertIn('human_g1k_v37_decoy.fasta', cmd)

	def test_clip_bam(self):
		clip_bam = ClipBam()
		clip_bam.input_bam = 'input.bam'
		clip_bam.output_bam = 'output.bam'
		clip_bam.reference_genome = 'human_g1k_v37_decoy.fasta'
		clip_bam.metrics_txt = 'metrics.txt'
		cmd = clip_bam.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('output.bam', cmd)
		self.assertIn('human_g1k_v37_decoy.fasta', cmd)
		self.assertIn('metrics.txt', cmd)






import unittest
from autoseq.tools.variantcalling import *
from autoseq.pipeline.clinseq import ClinseqPipeline


class TestVariantCalling(unittest.TestCase):
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
		self.test_clinseq_pipeline = ClinseqPipeline(self.sample_data, self.ref_data, {}, "/tmp/liqbio-test/",
											  "/media/clinseq/disk4/PROBIO/", 'FALSE')

	def test_haplotypecaller(self):
		haplotypecaller = HaplotypeCaller()
		haplotypecaller.input_bam = 'input.bam'
		haplotypecaller.reference_sequence = 'dummy_ref.fasta'
		haplotypecaller.dbSNP = 'dbsnp.vcf'
		haplotypecaller.interval_list = 'dummy_intervals'
		haplotypecaller.output = 'output.vcf'
		cmd = haplotypecaller.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('dbsnp.vcf', cmd)
		self.assertIn('dummy_intervals', cmd)
		self.assertIn('output.vcf', cmd)

	def test_strelka_germline(self):
		strelka_germline = StrelkaGermline()
		strelka_germline.input_bam = 'input.bam'
		strelka_germline.reference_sequence = 'dummy_ref.fasta'
		strelka_germline.target_bed = 'dummy_targets'
		strelka_germline.output_dir = 'output'
		strelka_germline.output_filtered_vcf = 'dummy_filtered.vcf'
		cmd = strelka_germline.command()
		self.assertIn('input.bam', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('dummy_targets', cmd)
		self.assertIn('output', cmd)
		self.assertIn('dummy_filtered.vcf', cmd)			

	def test_vardict(self):
		vardict = VarDict()
		vardict.input_tumor = "input_tumor.bam"
		vardict.input_normal = "input_normal.bam"
		vardict.reference_sequence = "dummy.fasta"
		vardict.reference_dict = {}
		vardict.target_bed = "dummy_targets.bed"
		vardict.tumorid = "tumorid"
		vardict.normalid = "normalid"
		vardict.output = "output.txt"
		cmd = vardict.command()
		self.assertIn('input_tumor.bam', cmd)
		self.assertIn('input_normal.bam', cmd)
		self.assertIn('dummy.fasta', cmd)
		self.assertIn('dummy_targets.bed', cmd)
		self.assertIn('output.txt', cmd)

	def test_strelka_somatic(self):
		strelka_somatic = StrelkaSomatic()
		strelka_somatic.input_tumor = "input_tumor.bam"
		strelka_somatic.input_normal = "input_normal.bam"
		strelka_somatic.input_indel_candidates = "dummy_indels.vcf"
		strelka_somatic.tumor_id = "tumor_id"
		strelka_somatic.normal_id = "normal_id"
		strelka_somatic.reference_sequence = "dummy_ref.fasta"
		strelka_somatic.target_bed = "dummy_targets.bed"
		strelka_somatic.output_dir = "output_dir"
		strelka_somatic.output_snvs_vcf = "dummy_snvs.vcf"
		strelka_somatic.output_indels_vcf = "dummy_indels.vcf"
		cmd = strelka_somatic.command()
		self.assertIn('input_tumor.bam', cmd)
		self.assertIn('input_normal.bam', cmd)
		self.assertIn('dummy_indels.vcf', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('dummy_targets.bed', cmd)
		self.assertIn('dummy_snvs.vcf', cmd)
		self.assertIn('dummy_indels.vcf', cmd)

	def test_mutect2_somatic(self):
		mutect = Mutect2Somatic()
		mutect.input_tumor = "input_tumor.bam"
		mutect.input_normal = "input_normal.bam"
		mutect.tumor_id = "tumor_id"
		mutect.normal_id = "normal_id"
		mutect.reference_sequence = "dummy_ref.fasta"
		mutect.bamout = "dummy_out.bam"
		mutect.output = "output.vcf"
		mutect.output_filtered = "output_filtered.vcf"
		mutect.interval_list = "dummy_intervals"
		cmd = mutect.command()
		self.assertIn('input_tumor.bam', cmd)
		self.assertIn('input_normal.bam', cmd)
		self.assertIn('tumor_id', cmd)
		self.assertIn('normal_id', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('dummy_out.bam', cmd)
		self.assertIn('dummy_intervals', cmd)
		self.assertIn('output.vcf', cmd)
		self.assertIn('output_filtered.vcf', cmd)

	def test_varscan(self):
		varscan = Varscan2Somatic()
		varscan.input_tumor = "input_tumor.bam"
		varscan.input_normal = "input_normal.bam"
		varscan.reference_sequence = "dummy_ref.fasta"
		varscan.normal_pileup = "normal_pileup"
		varscan.tumor_pileup = "tumor_pileup"
		varscan.output_snv = "output_snv.vcf"
		varscan.output_indel = "output_indel.vcf"
		varscan.output_somatic_snv = "output_somatic_snv.vcf"
		varscan.output_somatic_indel = "output_somatic_indel.vcf"
		cmd = varscan.command()
		self.assertIn('input_tumor.bam', cmd)
		self.assertIn('input_normal.bam', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('normal_pileup', cmd)
		self.assertIn('tumor_pileup', cmd)
		self.assertIn('output_indel.vcf', cmd)
		self.assertIn('output_snv.vcf', cmd)

	def test_vep(self):
		vep = VEP()
		vep.input_vcf = "input.vcf"
		vep.output_vcf = "output.vcf"
		vep.reference_sequence = "dummy.fasta"
		vep.vep_dir = "dummy_dir"
		vep.brca_exchange_vcf = "dummy_brca_exchange"
		cmd = vep.command()
		self.assertIn('input.vcf', cmd)
		self.assertIn('dummy.fasta', cmd)
		self.assertIn('dummy_dir', cmd)
		self.assertIn('output.vcf', cmd)
		self.assertIn('dummy_brca_exchange', cmd)

	def test_somaticseq(self):
		somaticseq = SomaticSeq()
		somaticseq.input_normal = "input_normal.bam"
		somaticseq.input_tumor = "input_tumor.bam"
		somaticseq.reference_sequence = "dummy_ref.fasta"
		somaticseq.input_mutect_vcf = "mutect.vcf"
		somaticseq.input_varscan_snv = "varscan_snv.vcf"
		somaticseq.input_varscan_indel = "varscan_indel.vcf"
		somaticseq.input_vardict_vcf = "vardict.vcf"
		somaticseq.input_strelka_snv = "strelka_snv.vcf"
		somaticseq.input_strelka_indel = "strelka_indel.vcf"
		somaticseq.output_dir = "output_dir"
		somaticseq.output_snv = "output_snv.vcf"
		somaticseq.output_indel = "output_indel.vcf"
		somaticseq.output_vcf = "output.vcf"
		cmd = somaticseq.command()
		self.assertIn('input_normal.bam', cmd)
		self.assertIn('input_tumor.bam', cmd)
		self.assertIn('dummy_ref.fasta', cmd)
		self.assertIn('mutect.vcf', cmd)
		self.assertIn('varscan_snv.vcf', cmd)
		self.assertIn('varscan_indel.vcf', cmd)
		self.assertIn('vardict.vcf', cmd)
		self.assertIn('strelka_snv.vcf', cmd)
		self.assertIn('strelka_indel.vcf', cmd)
		self.assertIn('output_dir', cmd)
		self.assertIn('output_snv.vcf', cmd)
		self.assertIn('output_indel.vcf', cmd)
		self.assertIn('output.vcf', cmd)

	def test_mergevcf(self):
		mergevcf = MergeVCF()
		mergevcf.input_vcf_hc = "input_hc.vcf"
		mergevcf.input_vcf_strelka = "input_strelka.vcf"
		mergevcf.output_vcf = "output.vcf"
		mergevcf.reference_genome = "dummy_ref.fasta"
		cmd = mergevcf.command()
		self.assertIn('input_hc.vcf', cmd)
		self.assertIn('input_strelka.vcf', cmd)
		self.assertIn('output.vcf', cmd)
		self.assertIn('dummy_ref.fasta', cmd)

	def test_generate_igvnav_input(self):
		generate_igvnav_input = GenerateIGVNavInput()
		generate_igvnav_input.input_vcf = "input.vcf"
		generate_igvnav_input.oncokb_db = "oncokb_db.txt"
		generate_igvnav_input.output = "output.txt"
		generate_igvnav_input.vcftype = "dummy_somatic"
		cmd = generate_igvnav_input.command()
		self.assertIn('input.vcf', cmd)
		self.assertIn('oncokb_db.txt', cmd)
		self.assertIn('output.txt', cmd)
		self.assertIn('dummy_somatic', cmd)

	def test_call_somatic_variants(self):
		num_jobs_before_call = len(self.test_clinseq_pipeline.graph.nodes())
		output_dict = call_somatic_variants(self.test_clinseq_pipeline, "test_cancer.bam", "test_normal.bam",
											self.test_tumor_capture, self.test_normal_capture, "progression",
											"test_outdir")
		self.assertListEqual(sorted(output_dict.keys()), sorted(['vardict','strelka_snvs', 'strelka_indels','mutect2', 'varscan_snv', 'varscan_indel']))
		num_jobs_after_call = len(self.test_clinseq_pipeline.graph.nodes())
		self.assertEquals(num_jobs_after_call, num_jobs_before_call + 4)

		






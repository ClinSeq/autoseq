from pypedream.job import required, Job, conditional, optional

__author__ = 'rebber'


class CreateContestVCFs(Job):
    """Runs create_contest_vcfs.py to generate a ContEst population allele frequency
    VCF input file, given an overall population VCF file and a pair of target region
    bed files."""
    
    def __init__(self):
        Job.__init__(self)
        self.input_target_regions_bed_1 = None
        self.input_target_regions_bed_2 = None
        self.input_population_vcf = None
        self.output = None
        self.jobname = "create_contest_vcfs"

    def command(self):
        return "create_contest_vcfs.py " + \
               required(" ", self.input_target_regions_bed_1) + \
               required(" ", self.input_target_regions_bed_2) + \
               required(" ", self.input_population_vcf) + \
               required("--output-filename ", self.output)


class ContEst(Job):
    """Runs ContEst to estimate contamination level in bam file "input_eval_bam"."""

    def __init__(self):
        Job.__init__(self)
        self.reference_genome = None
        self.input_eval_bam = None
        self.input_genotype_bam = None
        self.input_population_af_vcf = None
        self.output = None
        self.jobname = "contest"

    def command(self):
        min_genotype_ratio = "0.95"

        return "java -Xmx15g -jar /nfs/PROBIO/autoseq-scripts/GenomeAnalysisTK-3.5.jar -T ContEst " + \
            required("-R ", self.reference_genome) + \
            required("-I:eval ", self.input_eval_bam) + \
            required("-I:genotype ", self.input_genotype_bam) + \
            required("--popfile ", self.input_population_af_vcf) + \
            required("--min_genotype_ratio ", min_genotype_ratio) + \
            required("-o ", self.output)


class ContEstToContamCaveat(Job):
    """Runs script to convert ContEst output to JSON file with contamination QC
    estimatimate."""

    def __init__(self):
        Job.__init__(self)
        self.input_contest_results = None
        self.output = None
        self.jobname = "contest_to_contam_qc"

    def command(self):
        return "contest_to_contam_caveat.py " + \
            required(" ", self.input_contest_results) + \
            required("> ", self.output)
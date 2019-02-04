import os
from pypedream.job import Job, required, optional, conditional, repeat


class HeterzygoteConcordance(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.input_bam = None
        self.reference_sequence = None
        self.target_regions = None
        self.normalid = None
        self.filter_reads_with_N_cigar = True
        self.output = None
        self.jobname = "hzconcordance"

    def command(self):
        return "gatk-klevebring -T HeterozygoteConcordance " + \
               required("-R ", self.reference_sequence) + \
               required("-V ", self.input_vcf) + \
               required("-I ", self.input_bam) + \
               required("-sid ", self.normalid) + \
               optional("-L ", self.target_regions) + \
               conditional(self.filter_reads_with_N_cigar, "--filter_reads_with_N_cigar") + \
               required("-o ", self.output)


class FastQC(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.outdir = None
        self.extract = False
        self.jobname = "fastqc"

    def command(self):
        return "fastqc -o {} ".format(self.outdir) + \
               conditional(self.extract, "--extract") + \
               " --nogroup {}".format(self.input)


class MultiQC(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_files = None
        self.search_dir = None
        self.output = None
        self.report_title = None
        self.data_format = 'json'
        self.jobname = "multiqc"

    def command(self):
        required("input_files", self.input_files)
        required("output_base", self.output)
        required("dir_to_search", self.search_dir)

        basefn = os.path.basename(self.output)
        odir = os.path.dirname(self.output)
        return "multiqc " + \
               required("", self.search_dir) + \
               required("-o ", odir) + \
               optional("-n ", basefn) + \
               optional("-k ", self.data_format) + \
               optional("-i ", self.report_title) + \
               " --data-dir --zip-data-dir -v -f"


class SambambaDepth(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.targets_bed = None
        self.coverage_thresholds = [30, 50, 70, 100, 200, 300]
        self.output = None
        self.jobname = "sambamba-depth"

    def command(self):
        return "sambamba depth region" + \
               repeat("-T ", self.coverage_thresholds) + \
               required("--regions ", self.targets_bed) + \
               required("", self.input) + \
               required(">", self.output)


class BedtoolsCoverageHistogram(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bam = None
        self.input_bed = None
        self.tag = None
        self.output = None

    def command(self):
        tag_cmd = ""
        if self.tag:
            tag_cmd = "echo \"# {}\" >> {} \n".format(self.tag, self.output)
        return "echo \"# bedtools-coverage-hist: {}\"".format(self.input_bam) + \
               required("d>", self.output) + "\n" + \
               tag_cmd + \
               "bedtools coverage -hist " + \
               required("-a ", self.input_bed) + \
               required("-b ", self.input_bam) + \
               "|grep \"^all\" " + required(">> ", self.output)


class CoverageHistogram(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_bam = None
        self.input_bed = None
        self.min_basequal = None
        self.output = None

    def command(self):
        return "target_coverage_histogram.py " + \
               required("--targets ", self.input_bed) + \
               required(" ", self.input_bam) + \
               optional("--min_basequal ", self.min_basequal) + \
               required("> ", self.output)


class CoverageCaveat(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_histogram = None
        self.output = None
        self.high_thresh_fraction = 0.95
        self.high_thresh_fold_cov = 100
        self.low_thresh_fraction = 0.95
        self.low_thresh_fold_cov = 50

    def command(self):
        return "extract_coverage_caveat.py " + \
               required(" ", self.input_histogram) + \
               required("--high-thresh-fraction ", self.high_thresh_fraction) + \
               required("--high-thresh-fold-cov ", self.high_thresh_fold_cov) + \
               required("--low-thresh-fraction ", self.low_thresh_fraction) + \
               required("--low-thresh-fold-cov ", self.low_thresh_fold_cov) + \
               required("> ", self.output)

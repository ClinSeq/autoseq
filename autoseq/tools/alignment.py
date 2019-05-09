from pypedream.job import *
from pypedream.tools.unix import Cat

from autoseq.util.path import normpath
from autoseq.util.clinseq_barcode import *

__author__ = 'dankle'


class Bwa(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fastq1 = None  # input ports must start with "input"
        self.input_fastq2 = None
        self.input_reference_sequence = None
        self.mark_secondary = None
        self.remove_duplicates = True
        self.readgroup = None
        self.output = None  # output ports must start with "output", can be "output_metrics", "output", etc
        self.duplication_metrics = None
        self.jobname = "bwa"

    def command(self):
        bwalog = self.output + ".bwa.log"
        samblasterlog = self.output + ".samblaster.log"
        tmpprefix = "{}/{}".format(self.scratch, uuid.uuid4())

        return "bwa mem -M -v 1 " + \
               required("-R ", self.readgroup) + \
               optional("-t ", self.threads) + \
               required(" ", self.input_reference_sequence) + \
               required(" ", self.input_fastq1) + \
               optional("", self.input_fastq2) + \
               required("2>", bwalog) + \
               "| samblaster -M --addMateTags " + \
               conditional(self.remove_duplicates, "--removeDups") + \
               optional("--metricsFile ", self.duplication_metrics) + \
               required("2>", samblasterlog) + \
               "| samtools view -Sb -u - " + \
               "| samtools sort " + \
               required("-T ", tmpprefix) + \
               optional("-@ ", self.threads) + \
               required("-o ", self.output) + \
               " - " + \
               " && samtools index " + self.output + \
               " && cat {} {}".format(bwalog, samblasterlog) + \
               " && rm {} {}".format(bwalog, samblasterlog)


class Skewer(Job):
    def __init__(self):
        Job.__init__(self)
        self.input1 = None
        self.input2 = None
        self.output1 = None
        self.output2 = None
        self.stats = None
        self.jobname = "skewer"

    def command(self):
        if not self.output1.endswith(".gz") or not self.output2.endswith(".gz"):
            raise ValueError("Output files need to end with .gz")

        tmpdir = os.path.join(self.scratch, "skewer-" + str(uuid.uuid4()))
        prefix = "{}/skewer".format(tmpdir)

        if self.input2:
            out_fq1 = prefix + "-trimmed-pair1.fastq.gz"
            out_fq2 = prefix + "-trimmed-pair2.fastq.gz"
        else:
            out_fq1 = prefix + "-trimmed.fastq.gz"
            out_fq2 = ""

        out_stats = prefix + "-trimmed.log"

        mkdir_cmd = "mkdir -p {}".format(tmpdir)

        skewer_cmd = "skewer -z " + \
                     optional("-t ", self.threads) + " --quiet " + \
                     required("-o ", prefix) + \
                     required("", self.input1) + \
                     optional("", self.input2)
        copy_output_cmd = "cp " + out_fq1 + " " + self.output1 + \
            conditional(self.input2, " && cp " + out_fq2 + " " + self.output2)

        copy_stats_cmd = "cp " + out_stats + " " + self.stats
        rm_cmd = "rm -r {}".format(tmpdir)
        return " && ".join([mkdir_cmd, skewer_cmd, copy_output_cmd, copy_stats_cmd, rm_cmd])


class Realignment(Job):
    def __init__(self,):
        Job.__init__(self)
        self.input_bam = None
        self.output_bam = None
        self.reference_genome = None
        self.known_indel1 = None
        self.known_indel2 = None
        self.target_intervals = None
        self.target_region =  None
        self.jobname = "Realignment"

    def command(self):

        # creating target intervals for indel realignment 
        # Param: -L can be added to specify the genomic region
        target_creator_cmd = "java -jar /nfs/PROBIO/autoseq-scripts/GenomeAnalysisTK-3.5.jar " + \
                            " -T RealignerTargetCreator " + \
                            " -R " + self.reference_genome + \
                            " -known " + self.known_indel1 + \
			                " -allowPotentiallyMisencodedQuals " + \
                            " -L " + self.target_region + \
                            " -known " + self.known_indel2 + \
                            " -I " + self.input_bam + \
                            " -o " + self.target_intervals 

        realign_reads_cmd = "java -Xmx8G " + \
                            required("-Djava.io.tmpdir=", self.scratch) + \
                            " -jar /nfs/PROBIO/autoseq-scripts/GenomeAnalysisTK-3.5.jar " + \
                            " -T IndelRealigner " + \
                            " -R " + self.reference_genome + \
                            " -targetIntervals " + self.target_intervals + \
                            " -known " + self.known_indel1 + \
                            " -known " + self.known_indel2 + \
                            " -allowPotentiallyMisencodedQuals " + \
                            " -I " + self.input_bam + \
                            " -o " + self.output_bam

        return " && ".join([target_creator_cmd, realign_reads_cmd])


def align_library(pipeline, fq1_files, fq2_files, clinseq_barcode, ref, outdir, maxcores=1,
                  remove_duplicates=True):
    """
    Align fastq files for a PE library
    :param remove_duplicates:
    :param pipeline:
    :param fq1_files:
    :param fq2_files:
    :param lib:
    :param ref:
    :param outdir:
    :param maxcores:
    :return:
    """
    if not fq2_files:
        logging.debug("lib {} is SE".format(clinseq_barcode))
        return align_se(pipeline, fq1_files, clinseq_barcode, ref, outdir, maxcores, remove_duplicates)
    else:
        logging.debug("lib {} is PE".format(clinseq_barcode))
        return align_pe(pipeline, fq1_files, fq2_files, clinseq_barcode, ref, outdir, maxcores, remove_duplicates)


def align_se(pipeline, fq1_files, clinseq_barcode, ref, outdir, maxcores, remove_duplicates=True):
    """
    Align single end data
    :param pipeline:
    :param fq1_files:
    :param lib:
    :param ref:
    :param outdir:
    :param maxcores:
    :param remove_duplicates:
    :return:
    """
    logging.debug("Aligning files: {}".format(fq1_files))
    fq1_abs = [normpath(x) for x in fq1_files]
    fq1_trimmed = []
    for fq1 in fq1_abs:
        skewer = Skewer()
        skewer.input1 = fq1
        skewer.input2 = None
        skewer.output1 = outdir + "/skewer/{}".format(os.path.basename(fq1))
        skewer.output2 = outdir + "/skewer/unused-dummyfq2-{}".format(os.path.basename(fq1))
        skewer.stats = outdir + "/skewer/skewer-stats-{}.log".format(os.path.basename(fq1))
        skewer.threads = maxcores
        skewer.jobname = "skewer/{}".format(os.path.basename(fq1))
        skewer.scratch = pipeline.scratch
        skewer.is_intermediate = True
        fq1_trimmed.append(skewer.output1)
        pipeline.add(skewer)

    cat1 = Cat()
    cat1.input = fq1_trimmed
    cat1.output = outdir + "/skewer/{}_1.fastq.gz".format(clinseq_barcode)
    cat1.jobname = "cat/{}".format(clinseq_barcode)
    cat1.is_intermediate = False
    pipeline.add(cat1)

    bwa = Bwa()
    bwa.input_fastq1 = cat1.output
    bwa.input_reference_sequence = ref
    bwa.remove_duplicates = remove_duplicates

    library_id = parse_prep_id(clinseq_barcode)
    sample_string = compose_sample_str(extract_unique_capture(clinseq_barcode))

    bwa.readgroup = "\"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:ILLUMINA\"".format(\
        rg_id=clinseq_barcode, rg_sm=sample_string, rg_lb=library_id)

    bwa.threads = maxcores
    bwa.output = "{}/{}.bam".format(outdir, clinseq_barcode)
    bwa.scratch = pipeline.scratch
    bwa.jobname = "bwa/{}".format(clinseq_barcode)
    bwa.is_intermediate = False
    pipeline.add(bwa)

    return bwa.output


def align_pe(pipeline, fq1_files, fq2_files, clinseq_barcode, ref, outdir, maxcores=1, remove_duplicates=True):
    """
    align paired end data
    :param pipeline:
    :param fq1_files:
    :param fq2_files:
    :param lib:
    :param ref:
    :param outdir:
    :param maxcores:
    :param remove_duplicates:
    :return:
    """
    fq1_abs = [normpath(x) for x in fq1_files]
    fq2_abs = [normpath(x) for x in fq2_files]
    logging.debug("Trimming {} and {}".format(fq1_abs, fq2_abs))
    pairs = [(fq1_abs[k], fq2_abs[k]) for k in range(len(fq1_abs))]

    fq1_trimmed = []
    fq2_trimmed = []

    for fq1, fq2 in pairs:
        skewer = Skewer()
        skewer.input1 = fq1
        skewer.input2 = fq2
        skewer.output1 = outdir + "/skewer/libs/{}".format(os.path.basename(fq1))
        skewer.output2 = outdir + "/skewer/libs/{}".format(os.path.basename(fq2))
        skewer.stats = outdir + "/skewer/libs/skewer-stats-{}.log".format(os.path.basename(fq1))
        skewer.threads = maxcores
        skewer.jobname = "skewer/{}".format(os.path.basename(fq1))
        skewer.scratch = pipeline.scratch
        skewer.is_intermediate = True
        fq1_trimmed.append(skewer.output1)
        fq2_trimmed.append(skewer.output2)
        pipeline.add(skewer)

    cat1 = Cat()
    cat1.input = fq1_trimmed
    cat1.output = outdir + "/skewer/{}-concatenated_1.fastq.gz".format(clinseq_barcode)
    cat1.jobname = "cat1/{}".format(clinseq_barcode)
    cat1.is_intermediate = True
    pipeline.add(cat1)

    cat2 = Cat()
    cat2.input = fq2_trimmed
    cat2.jobname = "cat2/{}".format(clinseq_barcode)
    cat2.output = outdir + "/skewer/{}-concatenated_2.fastq.gz".format(clinseq_barcode)
    cat2.is_intermediate = True
    pipeline.add(cat2)

    bwa = Bwa()
    bwa.input_fastq1 = cat1.output
    bwa.input_fastq2 = cat2.output
    bwa.input_reference_sequence = ref
    bwa.remove_duplicates = remove_duplicates

    library_id = parse_prep_id(clinseq_barcode)
    sample_string = compose_sample_str(extract_unique_capture(clinseq_barcode))

    bwa.readgroup = "\"@RG\\tID:{rg_id}\\tSM:{rg_sm}\\tLB:{rg_lb}\\tPL:ILLUMINA\"".format(\
        rg_id=clinseq_barcode, rg_sm=sample_string, rg_lb=library_id)

    bwa.threads = maxcores
    bwa.output = "{}/{}.bam".format(outdir, clinseq_barcode)
    bwa.jobname = "bwa/{}".format(clinseq_barcode)
    bwa.scratch = pipeline.scratch
    bwa.is_intermediate = False
    pipeline.add(bwa)

    return bwa.output

def fq_trimming(pipeline, fq1_files, fq2_files, clinseq_barcode, ref, outdir, maxcores=1):
    fq1_abs = [normpath(x) for x in fq1_files]
    fq2_abs = [normpath(x) for x in fq2_files]
    logging.debug("Trimming {} and {}".format(fq1_abs, fq2_abs))
    pairs = [(fq1_abs[k], fq2_abs[k]) for k in range(len(fq1_abs))]

    fq1_trimmed = []
    fq2_trimmed = []

    for fq1, fq2 in pairs:
        skewer = Skewer()
        skewer.input1 = fq1
        skewer.input2 = fq2
        skewer.output1 = outdir + "/skewer/libs/{}".format(os.path.basename(fq1))
        skewer.output2 = outdir + "/skewer/libs/{}".format(os.path.basename(fq2))
        skewer.stats = outdir + "/skewer/libs/skewer-stats-{}.log".format(os.path.basename(fq1))
        skewer.threads = maxcores
        skewer.jobname = "skewer/{}".format(os.path.basename(fq1))
        skewer.scratch = pipeline.scratch
        skewer.is_intermediate = True
        fq1_trimmed.append(skewer.output1)
        fq2_trimmed.append(skewer.output2)
        pipeline.add(skewer)

    cat1 = Cat()
    cat1.input = fq1_trimmed
    cat1.output = outdir + "/skewer/{}-concatenated_1.fastq.gz".format(clinseq_barcode)
    cat1.jobname = "cat1/{}".format(clinseq_barcode)
    cat1.is_intermediate = True
    pipeline.add(cat1)

    cat2 = Cat()
    cat2.input = fq2_trimmed
    cat2.jobname = "cat2/{}".format(clinseq_barcode)
    cat2.output = outdir + "/skewer/{}-concatenated_2.fastq.gz".format(clinseq_barcode)
    cat2.is_intermediate = True
    pipeline.add(cat2)

    return cat1.output, cat2.output

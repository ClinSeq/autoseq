from pyjobmanager.job import Job, required


class BwaIndex(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fasta = None
        self.algorithm = "bwtsw"
        self.output = None
        self.jobname = "bwa-index"

    def command(self):
        return "bwa index" + \
               required("-a ", self.algorithm) + \
               required(" ", self.input_fasta)


class SamtoolsFaidx(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fasta = None
        self.output = None
        self.jobname = "samtools-faidx"

    def command(self):
        return "samtools faidx " + \
               required(" ", self.input_fasta)


class GenerateChrSizes(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fai = None
        self.output = None
        self.jobname = "make-chrsizes"

    def command(self):
        return "head -n 25 " + \
               required(" ", self.input_fai) + \
               "|cut -f 1-2 " + \
               required(">", self.output)

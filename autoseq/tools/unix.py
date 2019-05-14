from pypedream.job import Job, required, optional


# case class bwaIndex(ref:File) extends ExternalCommonArgs with  SingleCoreJob with OneDayJob {
#   @Input(doc="Input reference") val _ref: File = ref
#   @Output(doc="Input reference") val _out: File = ref + ".bwt"
#   def commandLine = bwaPath + " index -a bwtsw " + ref
#   this.jobName = "bwaIndex"
#   this.analysisName = this.jobName
#   this.isIntermediate = false
# }


class Copy(Job):
    def __init__(self, input_file, output_file):
        Job.__init__(self)
        self.input = input_file
        self.output = output_file
        self.jobname = "copy"

    def command(self):
        cmd = "cp " + \
               required(" ", self.input) + \
               required(" ", self.output)
        if self.input.endswith(".bam"):
            cmd += " && samtools index {}".format(self.output)
        
        if self.input.endswith(".vcf.gz"):
            cmd += " && tabix -p vcf {}".format(self.output)
        return cmd


class Gunzip(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "gunzip"

    def command(self):
        return "gzip -cd " + \
               required(" ", self.input) + \
               required(" > ", self.output)


class Curl(Job):
    def __init__(self):
        Job.__init__(self)
        self.remote = None
        self.output = None
        self.jobname = "curl"

    def command(self):
        return "curl " + \
               required(" ", self.remote) + \
               required(" > ", self.output)

class Bgzip(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.filetype = None
        self.jobname = "bgzip"

    def command(self):
        return "cat " + \
               required(" ", self.input) + \
               " | bgzip " + \
               required(" > ", self.output) + \
               " && tabix " + \
               optional("-p ", self.filetype) + \
               " {} ".format(self.output)

from pyjobmanager.job import required, Job


class FilterGTFChromosomes(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "filter-gtf-chrs"

    def command(self):
        return "gtf_filter_unused_chrs.py " + \
               required(" < ", self.input) + \
               required(" > ", self.output)


class FilterGTFGenes(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "filter-gtf-genes"

    def command(self):
        return """awk '{if($3=="gene"){print $0}}' """ + \
               required(" ", self.input) + \
               "| sort -k1,1 -k4,4n " + \
               required(" > ", self.output)


class GTF2GenePred(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "gtf2genepred"

    def command(self):
        return "gtfToGenePred -genePredExt " + required(" ", self.input) + \
               """ /dev/stdout |awk '{print $12"\t"$0}'|cut -f 1-11 """ + \
               required(" > ", self.output)

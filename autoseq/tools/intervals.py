from pyjobmanager.job import required, Job


class SlopIntervalList(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "slop-interval-list"

    def command(self):
        return "slopIntervalList.py " + \
               required(" < ", self.input) + \
               required(" > ", self.output)


class IntervalListToBed(Job):
    def __init__(self):
        Job.__init__(self)
        self.input = None
        self.output = None
        self.jobname = "interval-list-to-bed"

    def command(self):
        return "picard_interval_list_to_bed6_converter.py " + \
               required(" ", self.input) + \
               required(" ", self.output)

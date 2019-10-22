import uuid
from pyjobmanager.job import required, Job, conditional


class MsiSensor(Job):
    def __init__(self):
        Job.__init__(self)
        self.msi_sites = None
        self.input_normal_bam = None
        self.input_tumor_bam = None
        self.output = None
        self.jobname = "msisensor"

    def command(self):
        output_prefix = "{scratch}/msisensor-{uuid}".format(scratch=self.scratch,
                                                            uuid=uuid.uuid4())
        output_table = "{}".format(output_prefix)
        output_dis = "{}_dis".format(output_prefix)
        output_germline = "{}_germline".format(output_prefix)
        output_somatic = "{}_somatic".format(output_prefix)

        return "msisensor msi " + \
               required("-d ", self.msi_sites) + \
               required("-n ", self.input_normal_bam) + \
               required("-t ", self.input_tumor_bam) + \
               required("-o ", output_prefix) + \
               required("-b ", self.threads) + \
               " && cp {} {}".format(output_prefix, self.output) + \
               " && rm {} {} {} {}".format(output_table, output_dis,
                                           output_germline, output_somatic)


class MsiSensorScan(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fasta = None
        self.output = None
        self.homopolymers_only = True
        self.jobname = "msisensor-scan"

    def command(self):
        return "msisensor scan " + \
               required("-d ", self.input_fasta) + \
               required("-o ", self.output) + \
               conditional(self.homopolymers_only, " -p 1 ")


class IntersectMsiSites(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_msi_sites = None
        self.target_bed = None
        self.output_msi_sites = True
        self.jobname = "msi-intersect"

    def command(self):
        return "intersect-msi-sites.sh " + \
               required(" ", self.input_msi_sites) + \
               required(" ", self.target_bed) + \
               required(" ", self.output_msi_sites)


class Msings(Job):
    def __init__(self):
        Job.__init__(self)
        self.input_fasta = None
        self.outdir = None
        self.output = None
        self.msings_baseline = None
        self.msings_bed = None
        self.msings_intervals = None
        self.input_bam = None

        self.jobname = "msings"

    def command(self):
        return "run_msings.sh" + \
               required("-b ", self.msings_bed) + \
               required("-f ", self.input_fasta) + \
               required("-i ", self.msings_intervals) + \
               required("-n ", self.msings_baseline) + \
               required("-o ", self.outdir) + \
               required(" ", self.input_bam)

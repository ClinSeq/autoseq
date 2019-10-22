from pyjobmanager.job import required, Job
import uuid

__author__ = 'Thomas Whitington'


class MakeAllelicFractionTrack(Job):
    """
    Generate an IGV track representing variant allelic fraction information (bedGraph file).
    """

    def __init__(self):
        Job.__init__(self)
        self.input_vcf = None
        self.output_bedgraph = None
        self.jobname = "make_allelic_fraction_track"

    def command(self):
        return "generate_allelic_fraction_bedGraph.py " + \
               required("--output ", self.output_bedgraph) + \
               required(" ", self.input_vcf)


class MakeCNVkitTracks(Job):
    """
    Generate a IGV tracks representing the profile and segment information from a CNV-kit run.
    """

    def __init__(self):
        Job.__init__(self)
        self.input_cns = None
        self.input_cnr = None
        self.output_profile_bedgraph = None
        self.output_segments_bedgraph = None
        self.jobname = "make_cnvkit_tracks"

    def command(self):
        awk_cmd1 = "awk '$1 != \"chromosome\" {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}' %s > %s" % \
                   (self.input_cnr, self.output_profile_bedgraph)
        awk_cmd2 = "awk '$1 != \"chromosome\" {print $1\"\\t\"$2\"\\t\"$3\"\\t\"$5}' %s > %s" % \
                   (self.input_cns, self.output_segments_bedgraph)
        return "{} && {}".format(awk_cmd1, awk_cmd2)


class MakeQDNAseqTracks(Job):
    """
    Generate a IGV tracks representing the information from a QDNA-seq run.
    """

    def __init__(self):
        Job.__init__(self)
        self.input_qdnaseq_file = None
        self.output_segments_bedgraph = None
        self.output_copynumber_tdf = None
        self.output_readcount_tdf = None
        self.jobname = "make_qdnaseq_tracks"

    def command(self):
        copynumber_wig = "{scratch}/copynumber-{uuid}.wig".format(
            scratch=self.scratch, uuid=uuid.uuid4())
        readcount_wig = "{scratch}/readcount-{uuid}.wig".format(
            scratch=self.scratch, uuid=uuid.uuid4())

        qdnaseq_to_bedgraph_cmd = "qdnaseq_to_bedgraph.py {} {}".format(
            self.input_qdnaseq_file, self.output_segments_bedgraph)

        qdnaseq_to_wig_cmd = "qdnaseq_to_wig.py {} {} {}".format(
            self.input_qdnaseq_file, copynumber_wig, readcount_wig)
        igvtools_cmd1 = "igvtools toTDF {copynumber_wig} {copynumber_tdf} hg19".format(
            copynumber_wig=copynumber_wig, copynumber_tdf=self.output_copynumber_tdf
        )
        igvtools_cmd2 = "igvtools toTDF {readcount_wig} {readcount_tdf} hg19".format(
            readcount_wig=readcount_wig, readcount_tdf=self.output_readcount_tdf
        )

        return "{} && {} && {}".format(qdnaseq_to_bedgraph_cmd,
                                       qdnaseq_to_wig_cmd,
                                       igvtools_cmd1,
                                       igvtools_cmd2)
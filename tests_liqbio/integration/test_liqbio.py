import json
import os
import subprocess
import unittest

from genomicassertions.readassertions import ReadAssertions
from genomicassertions.variantassertions import VariantAssertions

#from autoseq.tests import liqbio_test_outdir


class TestLiqbio(unittest.TestCase, VariantAssertions, ReadAssertions):
    returncode = None
    tmpdir = None
    somatic_vcf = None
    outdir = "/media/clinseq/disk4/liqbio-test_data/test_output"
    libdir = "/media/clinseq/disk4/liqbio-test_data"
    ref_json = "/nfs/PROBIO/autoseq-genome/autoseq-genome.json"
    scratchdir = "/tmp/autoseq-integration-tests/liqbio"
    sample_json = "/media/clinseq/disk4/liqbio-test_data/sample.json"

    @classmethod
    def setUpClass(cls):
        cls.jobdb = os.path.join(cls.outdir, "jobdb.json")
        command = "autoseq " + \
            " --ref {} ". format(cls.ref_json) + \
            " --outdir {} ".format(cls.outdir) + \
            " --libdir {} ".format(cls.libdir) + \
            " --scratch {} --jobdb {} --cores 2 ".format(cls.scratchdir ,cls.jobdb) + \
            " liqbio {} ".format(cls.sample_json)
            
        import sys; print >> sys.stderr, command
        subprocess.check_call(command, shell=True)

    def test_jobdb(self):
        jobdb = json.load(open(self.jobdb, 'r'))
        self.assertEqual(set([job['status'] for job in jobdb['jobs']]),
                         {'COMPLETED'})

    def test_panel_bam_readgroups(self):
        bam = os.path.join(self.outdir, "bams", "CB", "LB-P-00150441-T-03289395-TP-CB-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'ID': 'LB-P-00150441-T-03289395-TP20180412-CB20180412',
                                                       'SM': 'LB-P-00150441-T-03289395',
                                                       'LB': 'TP20180412',
                                                       'PL': 'ILLUMINA'}])

        bam = os.path.join(self.outdir, "bams", "CB", "LB-P-00150441-N-03289093-TP-CB-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'LB': 'TP20180412',
                                                       'ID': 'LB-P-00150441-N-03289093-TP20180412-CB20180412',
                                                       'SM': 'LB-P-00150441-N-03289093-TP20180412',
                                                       'PL': 'ILLUMINA'}]
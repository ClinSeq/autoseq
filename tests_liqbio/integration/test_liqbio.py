import json
import os
import subprocess
import unittest

from genomicassertions.readassertions import ReadAssertions
from genomicassertions.variantassertions import VariantAssertions

from autoseq.tests import liqbio_test_outdir


class TestLiqbio(unittest.TestCase, VariantAssertions, ReadAssertions):
    returncode = None
    tmpdir = None
    somatic_vcf = None
    outdir = ""
    libdir = ""
    ref_json = "/tmp/test-genome/autoseq-genome.json"
    scratchdir = "/scratch/tmp/autoseq-integration-tests/liqbio"
    sample_json = "tests/liqbio-test-sample.json"

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

    def test_vardict_somatic(self):
        vcf = os.path.join(self.outdir, "variants",
                           "AL-P-NA12877-T-03098849-TD-TT-AL-P-NA12877-N-03098121-TD-TT.vardict-somatic.vcf.gz")

        # TP53 insertion: MU2185182, chr17:g.7578475->G
        # TP53 deletion: MU25947, chr17:g.7577558G>-
        # TP53 DNV: MU52971976, chr17:g.7574003GG>AA
        # PIK3CA hotspot E545K, MU5219, chr3:g.178936091G>A
        # PTEN hotspot R130Q, MU29098, chr10:g.89692905G>A
        # PTEN hotspot R233*, MU589331, chr10:g.89717672C>T
        # AR intron variant, MU50988553, chrX:g.66788924G>A

        self.assertVcfHasSample(vcf, 'AL-P-NA12877-N-03098121')  # N lib id is set
        self.assertVcfHasSample(vcf, 'AL-P-NA12877-T-03098849')  # T lib id is the merged library

        # deletion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7577557, 'AG', 'A')

        # insertion is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7578474, 'C', 'CG')

        # dnv (GG>AA) is called as a single event
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 17, 7574003, 'GG', 'AA')

        # PTEN hotspots are called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89692905, 'G', 'A')
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 10, 89717672, 'C', 'T')

        # PIK3CA hotspot is called
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, 3, 178936091, 'G', 'A')

    def test_panel_bam_readgroups(self):
        bam = os.path.join(self.outdir, "bams", "TT", "AL-P-NA12877-T-03098849-TD-TT-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'ID': 'AL-P-NA12877-T-03098849-TD1-TT1',
                                                       'SM': 'AL-P-NA12877-T-03098849',
                                                       'LB': 'TD1',
                                                       'PL': 'ILLUMINA'}])

        bam = os.path.join(self.outdir, "bams", "TT", "LB-P-NA12877-CFDNA-03098850-TD-TT-nodups.bam")
        self.assertBamHeaderElementEquals(bam, 'RG', [{'LB': 'TD1',
                                                       'ID': 'LB-P-NA12877-CFDNA-03098850-TD1-TT1',
                                                       'SM': 'LB-P-NA12877-CFDNA-03098850',
                                                       'PL': 'ILLUMINA'},
                                                      {'LB': 'TD1',
                                                       'ID': 'LB-P-NA12877-CFDNA-03098850-TD1-TT2',
                                                       'SM': 'LB-P-NA12877-CFDNA-03098850',
                                                       'PL': 'ILLUMINA'}]
                                          )

    def test_wgs_bam_coverage(self):
        bam = os.path.join(self.outdir, "bams/WG/",
                           "AL-P-NA12877-T-03098849-TD-WG-nodups.bam")
        self.assertBamHasCoverageAt(bam, 1, '3', 3617655)  # 1x coverage in that position
        self.assertBamHasCoverageAt(bam, 0, '3', 3618655)  # no coverage in that position

    def test_germline_vcf(self):
        vcf = os.path.join(self.outdir, "variants",
                           "AL-P-NA12877-N-03098121-TD-TT.freebayes-germline.vcf.gz")

        self.assertVcfHasSample(vcf, 'AL-P-NA12877-N-03098121')  # sample is the name of the normal lib id
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '3', 178925677, 'G', 'A')  # SNP
        self.assertVcfHasVariantWithChromPosRefAlt(vcf, '17', 7579643, 'CCCCCAGCCCTCCAGGT', 'C')  # deletion

    def test_germline_vcf_with_added_tumor_afs(self):
        vcf = os.path.join(self.outdir, "variants",
                           "AL-P-NA12877-N-03098121-TD-TT-and-AL-P-NA12877-T-03098849-TD-TT.germline-variants-with-somatic-afs.vcf.gz")

        self.assertVcfHasSample(vcf, 'AL-P-NA12877-N-03098121')
        self.assertVcfHasSample(vcf, 'AL-P-NA12877-T-03098849')
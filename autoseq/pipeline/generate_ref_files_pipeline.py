import json
import logging
import re

import sys
from pypedream.pipeline.pypedreampipeline import PypedreamPipeline
from pypedream.runners.shellrunner import Shellrunner

from autoseq.tools.genes import FilterGTFChromosomes, GTF2GenePred, FilterGTFGenes
from autoseq.tools.indexing import BwaIndex, SamtoolsFaidx, GenerateChrSizes
from autoseq.tools.intervals import SlopIntervalList, IntervalListToBed
from autoseq.tools.msi import MsiSensorScan, IntersectMsiSites
from autoseq.tools.picard import PicardCreateSequenceDictionary
from autoseq.tools.qc import *
from autoseq.tools.unix import Gunzip, Curl, Copy, Bgzip
from autoseq.tools.variantcalling import VcfFilter, CurlSplitAndLeftAlign, InstallVep
from autoseq.util.path import stripsuffix, normpath

__author__ = 'dankle'


class GenerateRefFilesPipeline(PypedreamPipeline):
    outdir = None
    maxcores = None

    def __init__(self, genome_resources, outdir, maxcores=1, runner=Shellrunner()):
        PypedreamPipeline.__init__(self, normpath(outdir), runner=runner)

        self.genome_resources = genome_resources
        self.input_reference_sequence = "{}/human_g1k_v37_decoy.fasta.gz".format(genome_resources)
        self.cosmic_vcf = "{}/CosmicCodingMuts_v71.vcf.gz".format(genome_resources)
        self.qdnaseq_background = "{}/qdnaseq_background.Rdata".format(genome_resources)
        self.swegene_common_vcf = "{}/swegen_common.vcf.gz".format(genome_resources)
        self.thousand_genome_vcf = "{}/1000G_phase1.indels.b37.vcf.gz".format(genome_resources)
        self.mills_and_1000g_gold_standard = "{}/Mills_and_1000G_gold_standard.indels.b37.vcf.gz".format(genome_resources)
        self.brca_exchange = "{}/BrcaExchangeClinvar_15Jan2019_v26_hg19.vcf.gz".format(genome_resources)
        self.oncokb = "{}/OncoKB_6Mar19_v1.9.txt".format(genome_resources)
        self.outdir = outdir
        self.maxcores = maxcores
        self.reference_data = dict()

        self.exac_remote = "ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz"
        self.dbsnp_remote = "ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b149_GRCh37p13/VCF/All_20161121.vcf.gz"
        self.clinvar_remote = "ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/archive_1.0/2016/clinvar_20160203.vcf.gz"
        self.icgc_somatic_remote = "https://dcc.icgc.org/api/v1/download?fn=/release_20/Summary/simple_somatic_mutation.aggregated.vcf.gz"
        self.ensembl_version = "75"
        self.ensembl_gtf_remote = "ftp://ftp.ensembl.org/pub/release-" + self.ensembl_version + \
                                  "/gtf/homo_sapiens/Homo_sapiens.GRCh37." + self.ensembl_version + ".gtf.gz"
        self.mitranscriptome_remote = "http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz"

        self.prepare_reference_genome()
        self.prepare_genes()
        self.prepare_sveffect_regions()
        self.prepare_intervals()
        self.prepare_variants()

        fetch_vep_cache = InstallVep()
        fetch_vep_cache.output_dir = "{}/vep/".format(self.outdir)
        #self.add(fetch_vep_cache)

        self.reference_data['vep_dir'] = fetch_vep_cache.output_dir

        self.make_ref_paths_relative()

        with open("{}/autoseq-genome.json".format(self.outdir), "w") as output_file:
            json.dump(self.reference_data, output_file, indent=4, sort_keys=True)

    def prepare_sveffect_regions(self):
        for regions_name in ["ar_regions", "ts_regions", "fusion_regions"]:
            file_full_path = "{}/{}.bed".format(
                self.genome_resources,
                regions_name,
            )

            copy_regions = Copy(input_file=file_full_path,
                                output_file="{}/intervals/{}".format(
                                    self.outdir,
                                    os.path.basename(file_full_path),
                                ))

            self.reference_data[regions_name] = copy_regions.output

            self.add(copy_regions)


    def prepare_variants(self):
        curl_dbsnp = CurlSplitAndLeftAlign()
        curl_dbsnp.input_reference_sequence = self.reference_data['reference_genome']
        curl_dbsnp.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_dbsnp.remote = self.dbsnp_remote
        curl_dbsnp.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.dbsnp_remote))
        curl_dbsnp.is_intermediate = True
        self.add(curl_dbsnp)

        filter_dbsnp = VcfFilter()
        filter_dbsnp.input = curl_dbsnp.output
        filter_dbsnp.filter = "\"! ( SAO = 3 | SAO = 2 )\""
        filter_dbsnp.output = "{}/variants/dbsnp142-germline-only.vcf.gz".format(self.outdir)
        self.add(filter_dbsnp)

        curl_cosmic = CurlSplitAndLeftAlign()
        curl_cosmic.input_reference_sequence = self.reference_data['reference_genome']
        curl_cosmic.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_cosmic.remote = "file://" + self.cosmic_vcf
        curl_cosmic.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.cosmic_vcf))
        self.add(curl_cosmic)

        curl_clinvar = CurlSplitAndLeftAlign()
        curl_clinvar.input_reference_sequence = self.reference_data['reference_genome']
        curl_clinvar.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_clinvar.remote = self.clinvar_remote
        curl_clinvar.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.clinvar_remote))
        self.add(curl_clinvar)

        curl_exac = CurlSplitAndLeftAlign()
        curl_exac.input_reference_sequence = self.reference_data['reference_genome']
        curl_exac.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_exac.remote = self.exac_remote
        curl_exac.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.exac_remote))
        self.add(curl_exac)

        curl_icgc = CurlSplitAndLeftAlign()
        curl_icgc.input_reference_sequence = self.reference_data['reference_genome']
        curl_icgc.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_icgc.remote = self.icgc_somatic_remote
        curl_icgc.output = "{}/variants/{}".format(self.outdir,
                                                   "icgc_release_20_simple_somatic_mutation.aggregated.vcf.gz")
        self.add(curl_icgc)

        curl_swegene = CurlSplitAndLeftAlign()
        curl_swegene.input_reference_sequence = self.reference_data['reference_genome']
        curl_swegene.input_reference_sequence_fai = self.reference_data['reference_genome'] + ".fai"
        curl_swegene.remote = "file://" + self.swegene_common_vcf
        curl_swegene.output = "{}/variants/{}".format(self.outdir, os.path.basename(self.swegene_common_vcf))
        self.add(curl_swegene)

        copy_thousand_genome = Copy(input_file=self.thousand_genome_vcf,
                                    output_file="{}/variants/{}".format(self.outdir, os.path.basename(self.thousand_genome_vcf)))
        self.add(copy_thousand_genome)

        copy_mills_and_1000g = Copy(input_file=self.mills_and_1000g_gold_standard,
                                    output_file="{}/variants/{}".format(self.outdir, os.path.basename(self.mills_and_1000g_gold_standard)))
        self.add(copy_mills_and_1000g)

        copy_brca_exchange = Copy(input_file=self.brca_exchange,
                                output_file="{}/variants/{}".format(self.outdir, os.path.basename(self.brca_exchange)))
        self.add(copy_brca_exchange)

        copy_oncokb = Copy(input_file=self.oncokb,
                        output_file="{}/variants/{}".format(self.outdir, os.path.basename(self.oncokb)))
        self.add(copy_oncokb)

        self.reference_data['dbSNP'] = filter_dbsnp.output
        self.reference_data['cosmic'] = curl_cosmic.output
        self.reference_data['exac'] = curl_exac.output
        self.reference_data['clinvar'] = curl_clinvar.output
        self.reference_data['icgc'] = curl_icgc.output
        self.reference_data['swegene_common'] = curl_swegene.output
        self.reference_data['1KG'] = copy_thousand_genome.output
        self.reference_data['Mills_and_1KG_gold_standard'] = copy_mills_and_1000g.output
        self.reference_data['brca_exchange'] = copy_brca_exchange.output
        self.reference_data['oncokb'] = copy_oncokb.output
        

    def prepare_cnvkit(self, cnv_kit_ref_filename):
        """

        :param cnv_kit_ref_filename: String name of a cnvkit reference file (either *.cnn
        or *.cnvkitref.txt), to be registered in self.ref_data.
        """

        file_full_path = "{}/target_intervals/{}".format(self.genome_resources, cnv_kit_ref_filename)

        # Extract the capture+library+sampletype strings:
        capture_library_sampletype = cnv_kit_ref_filename.split(".")[:3]

        copy_cnvkit_ref = Copy(input_file=file_full_path,
                               output_file="{}/intervals/targets/{}".format(self.outdir,
                                                                            os.path.basename(file_full_path))
                              )
        self.add(copy_cnvkit_ref)

        capture_name = capture_library_sampletype[0]
        library_kit_name = capture_library_sampletype[1]
        sample_type = capture_library_sampletype[2]

        ref_type = "cnvkit-fix"
        if cnv_kit_ref_filename.endswith(("cnn")):
            ref_type = "cnvkit-ref"

        # FIXME: Ugly; refactor. This registers the cnvkit reference file copy in the autoseq genome dictionary:
        if ref_type not in self.reference_data['targets'][capture_name]:
            self.reference_data['targets'][capture_name][ref_type] = {}
        if library_kit_name not in self.reference_data['targets'][capture_name][ref_type]:
            self.reference_data['targets'][capture_name][ref_type][library_kit_name] = {}
        self.reference_data['targets'][capture_name][ref_type][library_kit_name][sample_type] = \
            copy_cnvkit_ref.output

    def prepare_msings(self, filename_base, capture_name):
        """
        Setup the copying of the relevant msings parameter files, for the given base
        filename.

        :param filename_base: String suffix for all potential msings parameter files.
        :param capture_name: String capture name. Note: this information is also contained
            in the filename_base parameter; should be refactored.
        """

        for msings_extn in ["baseline", "bed", "msi_intervals"]:
            msings_ref_file = filename_base + ".msings." + msings_extn
            if os.path.exists(msings_ref_file):
                copy_msings_ref = Copy(input_file=msings_ref_file,
                                       output_file="{}/intervals/targets/{}".format(self.outdir,
                                                                                    os.path.basename(
                                                                                        msings_ref_file))
                                       )
                self.add(copy_msings_ref)
                self.reference_data['targets'][capture_name]['msings-' + msings_extn] = copy_msings_ref.output
            else:
                self.reference_data['targets'][capture_name]['msings-' + msings_extn] = None

    # FIXME: The prepare_intervals() method is becoming very unwieldy. Consider refactoring.
    def prepare_intervals(self):
        self.reference_data['targets'] = {}
        target_intervals_dir = "{}/target_intervals/".format(self.genome_resources)
        input_files = [f for f in os.listdir(target_intervals_dir) if f.endswith(".interval_list")]

        scan_for_microsatellites = MsiSensorScan()
        scan_for_microsatellites.input_fasta = self.reference_data['reference_genome']
        scan_for_microsatellites.homopolymers_only = True
        scan_for_microsatellites.output = "{}/intervals/msisensor-microsatellites.tsv".format(self.outdir)
        self.add(scan_for_microsatellites)

        for f in input_files:
            file_full_path = "{}/target_intervals/{}".format(self.genome_resources, f)
            logging.debug("Parsing intervals file {}".format(file_full_path))
            capture_name = stripsuffix(f, ".interval_list")
            self.reference_data['targets'][capture_name] = {}

            copy_file = Copy(input_file=file_full_path,
                             output_file="{}/intervals/targets/{}".format(self.outdir,
                                                                          os.path.basename(file_full_path)))
            self.add(copy_file)

            slop_interval_list = SlopIntervalList()
            slop_interval_list.input = copy_file.output
            slop_interval_list.output = stripsuffix(copy_file.output, ".interval_list") + ".slopped20.interval_list"
            self.add(slop_interval_list)

            interval_list_to_bed = IntervalListToBed()
            interval_list_to_bed.input = slop_interval_list.output
            interval_list_to_bed.output = stripsuffix(slop_interval_list.output, ".interval_list") + ".bed"
            self.add(interval_list_to_bed)
            
            gzip_bed = Bgzip()
            gzip_bed.input = interval_list_to_bed.output
            gzip_bed.output = interval_list_to_bed.output + ".gz"
            gzip_bed.filetype = "bed"
            self.add(gzip_bed)


            intersect_msi = IntersectMsiSites()
            intersect_msi.input_msi_sites = scan_for_microsatellites.output
            intersect_msi.target_bed = interval_list_to_bed.output
            intersect_msi.output_msi_sites = stripsuffix(interval_list_to_bed.output, ".bed") + ".msisites.tsv"
            self.add(intersect_msi)

            self.prepare_msings(stripsuffix(file_full_path, ".interval_list"), capture_name)

            self.reference_data['targets'][capture_name]['blacklist-bed'] = None
            blacklist_bed = stripsuffix(file_full_path, ".interval_list") + ".blacklist.bed"
            if os.path.exists(blacklist_bed):
                blacklist_copy = Copy(input_file=blacklist_bed,
                                      output_file="{}/intervals/targets/{}".format(
                                          self.outdir,
                                          os.path.basename(blacklist_bed),
                                      ))
                self.add(blacklist_copy)
                self.reference_data['targets'][capture_name]['blacklist-bed'] = blacklist_copy.output

            purecn_targets_file = stripsuffix(file_full_path, ".interval_list") + ".purecn.txt"
            if os.path.exists(purecn_targets_file):
                copy_purecn_targets = Copy(input_file=purecn_targets_file,
                                       output_file="{}/intervals/targets/{}".format(self.outdir,
                                                                                    os.path.basename(
                                                                                        purecn_targets_file))
                                       )
                self.add(copy_purecn_targets)
                self.reference_data['targets'][capture_name]['purecn_targets'] = copy_purecn_targets.output
            else:
                self.reference_data['targets'][capture_name]['purecn_targets'] = None

            self.reference_data['targets'][capture_name]['targets-interval_list'] = copy_file.output
            self.reference_data['targets'][capture_name]['targets-interval_list-slopped20'] = slop_interval_list.output
            self.reference_data['targets'][capture_name]['targets-bed-slopped20'] = interval_list_to_bed.output
            self.reference_data['targets'][capture_name]['targets-bed-slopped20-gz'] = gzip_bed.output
            self.reference_data['targets'][capture_name]['msisites'] = intersect_msi.output_msi_sites

        # Find all .cnn files and copy + register them for use in cnv kit:
        for f in [f for f in os.listdir(target_intervals_dir) if (f.endswith(".cnn") or "cnvkit-fix" in f)]:
            self.prepare_cnvkit(f)


    def prepare_genes(self):
        curl_ensembl_gtf = Curl()
        curl_ensembl_gtf.remote = self.ensembl_gtf_remote
        curl_ensembl_gtf.output = "{}/genes/{}".format(self.outdir, os.path.basename(self.ensembl_gtf_remote))
        curl_ensembl_gtf.jobname = "curl-ensembl-gtf"
        curl_ensembl_gtf.is_intermediate = True
        self.add(curl_ensembl_gtf)

        gunzip_ensembl_gtf = Gunzip()
        gunzip_ensembl_gtf.input = curl_ensembl_gtf.output
        gunzip_ensembl_gtf.output = stripsuffix(curl_ensembl_gtf.output, ".gz")
        gunzip_ensembl_gtf.is_intermediate = True
        self.add(gunzip_ensembl_gtf)

        filt_ensembl_gtf_chrs = FilterGTFChromosomes()
        filt_ensembl_gtf_chrs.input = gunzip_ensembl_gtf.output
        filt_ensembl_gtf_chrs.output = stripsuffix(gunzip_ensembl_gtf.output, ".gtf") + ".filtered.gtf"
        self.add(filt_ensembl_gtf_chrs)

        gtf2genepred_ensembl = GTF2GenePred()
        gtf2genepred_ensembl.input = filt_ensembl_gtf_chrs.output
        gtf2genepred_ensembl.output = stripsuffix(filt_ensembl_gtf_chrs.output, ".gtf") + ".genepred"
        self.add(gtf2genepred_ensembl)

        filt_genes_ensembl_gtf_genes = FilterGTFGenes()
        filt_genes_ensembl_gtf_genes.input = filt_ensembl_gtf_chrs.output
        filt_genes_ensembl_gtf_genes.output = stripsuffix(filt_ensembl_gtf_chrs.output, ".gtf") + ".genes-only.gtf"
        self.add(filt_genes_ensembl_gtf_genes)

        self.reference_data['ensemblVersion'] = self.ensembl_version
        self.reference_data['genesGtf'] = filt_ensembl_gtf_chrs.output
        self.reference_data['genesGenePred'] = gtf2genepred_ensembl.output
        self.reference_data['genesGtfGenesOnly'] = filt_genes_ensembl_gtf_genes.output

    def prepare_reference_genome(self):
        genome_unzipped = stripsuffix(os.path.basename(self.input_reference_sequence), ".gz")

        gunzip_ref = Gunzip()
        gunzip_ref.input = self.input_reference_sequence
        gunzip_ref.output = "{}/genome/{}".format(self.outdir, genome_unzipped)
        self.add(gunzip_ref)

        copy_ref_to_bwa = Copy(input_file=gunzip_ref.output,
                               output_file="{}/bwa/{}".format(self.outdir, os.path.basename(gunzip_ref.output)))
        self.add(copy_ref_to_bwa)

        bwa_index = BwaIndex()
        bwa_index.input_fasta = copy_ref_to_bwa.output
        bwa_index.output = copy_ref_to_bwa.output + ".bwt"
        bwa_index.algorithm = "bwtsw"
        self.add(bwa_index)

        create_dict = PicardCreateSequenceDictionary()
        create_dict.input = gunzip_ref.output
        create_dict.output_dict = gunzip_ref.output.replace(".fasta", "") + ".dict"
        self.add(create_dict)

        samtools_faidx = SamtoolsFaidx()
        samtools_faidx.input_fasta = gunzip_ref.output
        samtools_faidx.output = gunzip_ref.output + ".fai"
        self.add(samtools_faidx)

        create_chrsizes = GenerateChrSizes()
        create_chrsizes.input_fai = samtools_faidx.output
        create_chrsizes.output = gunzip_ref.output.replace(".fasta", "") + ".chrsizes.txt"
        self.add(create_chrsizes)

        copy_qdnaseq_bg = Copy(input_file=self.qdnaseq_background,
                               output_file="{}/genome/{}".format(self.outdir,
                                                                 os.path.basename(self.qdnaseq_background)))
        self.add(copy_qdnaseq_bg)

        self.reference_data['reference_genome'] = gunzip_ref.output
        self.reference_data['reference_dict'] = create_dict.output_dict
        self.reference_data['chrsizes'] = create_chrsizes.output
        self.reference_data['bwaIndex'] = bwa_index.input_fasta
        self.reference_data['qdnaseq_background'] = copy_qdnaseq_bg.output

    def make_ref_paths_relative(self):
        """Recursively traverse a given dictionary and make paths relative"""

        def make_paths_relative(d):
            for k, v in d.items():
                if isinstance(v, dict):
                    make_paths_relative(v)
                else:
                    if v and self.outdir in v:
                        d[k] = os.path.relpath(v, self.outdir)
                        
        make_paths_relative(self.reference_data)


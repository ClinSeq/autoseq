# autoseq

[![Coverage Status](https://coveralls.io/repos/github/ClinSeq/autoseq/badge.svg?branch=master)](https://coveralls.io/github/ClinSeq/autoseq?branch=master)

## Clinseq barcodes

Each sample+preparation+capture item should have a corresponding barcode with the format `PROJECT-SDID-TYPE-SAMPLEID-PREPID-CAPTUREID` where:

* `PROJECT` is a two-letter short project designator. One of `AL` (alascca), `LB` (liquid biopspy) and `OT` (other)
* `SDID` is an identifier for a single individual. It must match the pattern `P-[a-zA-Z0-9]+` (*NOTE:* This necessitates an additional "-" within this field).
* `TYPE` is the sample type, one of `T` (tumor), `N` (normal) and `CFDNA` (ctDNA)
* `SAMPLEID` identifies a single biological sample, for example piece of a tumor or a single tube of plasma. It must match the pattern `[a-zA-Z0-9]+`.
* `PREPID` specifies the library preparation kit used. It must be a two-letter shortname followed by a string matching `[0-9]+`, which can be used to indicate the date on which the prep was performed. The date string should *preferably* be in the format `YYYYMMDDHHMM`. For example, `201701241540` would indicate year 2017, January 24th, at 15:40.
* `CAPTUREID` specifies the capture that was performed on the library (if any). It must match either `WGS` (indicating that no capture was performed), or else a two-letter shortname indicating the capture kit used, followed by a string matching `[0-9]+`, which can be used to indicate the date on which the capture was performed. The date should *preferably* be in the format `YYYYMMDDHHMM`.

*NOTE:* The combination `SDID-TYPE-SAMPLEID` must uniquely identify a single sample.

*NOTE:* A clinseq barcode is not garuanteed to uniquely specify a single sample+library+capture item, but in practice it should be unique if precise preparation and capture times are included within the `PREPID` and `CAPTUREID` fields.

### Allowed Prep IDs

Autoseq know about the following preparation methods: 

* `BN` = `BIOO_NEXTFLEX`
* `KH` = `KAPA_HYPERPREP`
* `KP` = `KAPA_HYPERPLUS`
* `TD` = `THRUPLEX_DNASEQ`
* `TP` = `THRUPLEX_PLASMASEQ`
* `TF` = `THRUPLEX_FD`
* `TS` = `TRUSEQ_RNA`
* `NN` = `NEBNEXT_RNA`
* `VI` = `VILO_RNA`            

### Allowed Capture IDs

Autoseq knows about the following capture kits:

* `CS` = `clinseq_v3_targets`
* `CZ` = `clinseq_v4`
* `EX` = `EXOMEV3`
* `EO` = `EXOMEV1`
* `RF` = `fusion_v1`
* `CC` = `core_design`
* `CD` = `discovery_coho`
* `CB` = `big_design`
* `TT` = `test-regions`
* `CM` = `monitor`
* `CP` = `progression`

## Runners

Autoseq can use any of the runners implemented in pypedream, `shellrunner` (default), `localqrunner` or `slurmrunner`.

# General options

`--libdir` is the directory where the libraries live. Each library should have its own subdirectory where fastq.gz files can be placed. Autoseq recoginzes files on the format `_1.fastq.gz/_2.fastq.gz`.


# LiqBio pipeline

The Liquid Biopsy pipeline is invoked by

~~~
autoseq --ref ref.json --outdir /path/to/outdir --jobdb jobdb.json --cores 5 --runner_name slurmrunner --libdir /path/to/libdir liqbio sample.json
~~~

The `sample.json` file has the format

~~~json
{
    "sdid": "NA12877",
    "panel": {
        "T": "NA12877-T-03098849-TD1-TT1",
        "N": "NA12877-N-03098121-TD1-TT1",
        "CFDNA": ["NA12877-CFDNA-03098850-TD1-TT1", "NA12877-CFDNA-03098850-TD2-TT1"]
    },
    "wgs": {
        "T": "NA12877-T-03098849-TD1-WGS",
        "N": "NA12877-N-03098121-TD1-WGS",
        "CFDNA": ["NA12877-CFDNA-03098850-TD1-WGS"]
    }
}

~~~

In this file, a single tumor and normal sample is allowed, but multiple plasma samples. If no tumor or normal sample is avaialble, they can be set to `null`, but if no plasma samples are available, it should be set to `[]` (empty list), for example `"CFDNA": []`.

For the plasma samples, merging of libraries will take place before calling. On alignment, the `@RG` tag will be set as follows: 

* `ID` = `SDID-TYPE-SAMPLEID-PREPID-CAPTUREID`
* `LB` = `SDID-TYPE-SAMPLEID-PREPID`
* `SM` = `SDID-TYPE-SAMPLEID`

Of note is that the library tag (`LB`) does not include the `CAPTUREID` part, to ensure that PCR duplicates are removed correctly. 

If a single prepared samples is exposed to capture twice, to create the libraries `NA12877-T-49-TD1-TT1` and `NA12877-T-49-TD1-TT2` (note different digits in the capture id), read pairs being identical between the two libraries should be considered duplicates since the sample was split after the final PCR step. Therefore, the `LB` for these libraries is set to `NA12877-T-49-TD1`. After merging the bam files, removal of PCR duplicates is done using Picard MarkDuplicates, which will do the right thing. 

# Automated testing on travis-ci

For automated testing, a test reference genome and a test datas set with relevant data are supplied. 

### Reference genome

The test reference genome and assets is available for download at `https://export.uppmax.uu.se/b2010040/test-genome.tar.gz`. This archive contains a sliced version of a full set of genome files for autoseq, including various key genes. 

The whole chromosomes 3, 10, 17, X and Y are selected, after which everything except the following regions have been masked (to speed up alignment):

~~~
3	178863388	179014224	PIK3CA_150k
10	83068546	96283182	PTEN_13M
17	7558477	7589399	TP53_30k
X	66782057	66796840	14k_AR_exon
Y	6810425	6825985	15k_on_Y
~~~

From these regions, key exons and various other regions have been selected to mimic a small exome. 

### The Test Dataset

A sythetic tumor/normal/plasma dataset has been created for testing purpuses. From the illumina platinum 200x WGS sample from NA12877, read pairs from the seleted targets have been extracted. These reads have then been randomly assigned to create a virtual normal sample with ≈50x coverage, and remaining reads (≈150x coverage) have been put aside. To create a virtual tumor and a virtual plasma sample, variants have been spiked into the 150x data in the following positions: 

* TP53 insertion: MU2185182, chr17:g.7578475->G
* TP53 deletion: MU25947, chr17:g.7577558G>-
* TP53 DNV: MU52971976, chr17:g.7574003GG>AA
* PIK3CA hotspot E545K, MU5219, chr3:g.178936091G>A
* PTEN hotspot R130Q, MU29098, chr10:g.89692905G>A
* PTEN hotspot R233*, MU589331, chr10:g.89717672C>T
* AR intron variant, MU50988553, chrX:g.66788924G>A

In the virtual tumor, the target variant allele fraction (VAF) is 30% and in the virtual plasma sample the target VAF is 20%. 

The variants have been selected from ICGC simple somatic mutations v20 with the aim to cover common small variants, including SNVs, deletions, insertions and DNVs. Note that the tests does not address the issue of global sensitivity and PPV of the pipeline, but are only intented to ensure that variants of all kinds are detected by the pipeline. 

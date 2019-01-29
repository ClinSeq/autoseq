import json
import logging
import os
import sys
import time

import click

from autoseq.pipeline.liqbio import LiqBioPipeline
from autoseq.util.clinseq_barcode import extract_clinseq_barcodes, convert_barcodes_to_sampledict, validate_clinseq_barcodes
from autoseq.util.path import mkdir


@click.command()
@click.argument('sample', type=click.File('r'))
@click.pass_context
def liqbio(ctx, sample):
    logging.info("Running Liquid Biopsy pipeline")
    logging.info("Sample is {}".format(sample))

    logging.debug("Reading sample config from {}".format(sample))
    sampledata = json.load(sample)

    if ctx.obj['jobdb']:
        mkdir(os.path.dirname(ctx.obj['jobdb']))
    
    ctx.obj['pipeline'] = LiqBioPipeline(sampledata=sampledata,
                                         refdata=ctx.obj['refdata'],
                                         job_params=ctx.obj['job_params'],
                                         outdir=ctx.obj['outdir'],
                                         libdir=ctx.obj['libdir'],
                                         maxcores=ctx.obj['cores'],
                                         runner=ctx.obj['runner'],
                                         jobdb=ctx.obj['jobdb'],
                                         dot_file=ctx.obj['dot_file'],
                                         umi=ctx.obj['umi'],
                                         scratch=ctx.obj['scratch'])

    # start main analysis
    ctx.obj['pipeline'].start()
    #
    logging.info("Waiting for pipeline to finish.")
    while ctx.obj['pipeline'].is_alive():
        logging.debug("Waiting for LiqBioPipeline")
        time.sleep(5)

    # # return_code from run_pipeline() will be != 0 if the pipeline fails
    sys.exit(ctx.obj['pipeline'].exitcode)


@click.command()
@click.option('--outdir', required=True, help="directory to write config files")
@click.argument('barcodes-filename', type=str)
@click.pass_context
def liqbio_prepare(ctx, outdir, barcodes_filename):
    logging.info("Extracting clinseq barcodes from input file: " + barcodes_filename)
    clinseq_barcodes = extract_clinseq_barcodes(barcodes_filename)

    logging.info("Validating all clinseq barcodes.")
    validate_clinseq_barcodes(clinseq_barcodes)

    logging.info("Generating sample dictionary from the input clinseq barcode strings.")
    sample_dict = convert_barcodes_to_sampledict(clinseq_barcodes)
    for sdid in sample_dict:
        fn = "{}/{}.json".format(outdir, sdid)
        with open(fn, 'w') as f:
            json.dump(sample_dict[sdid], f, sort_keys=True, indent=4)

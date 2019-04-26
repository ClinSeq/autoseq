import json
import logging
import os
import signal

import click
from pypedream import runners

from .alascca import alascca as alascca_cmd
from .liqbio import liqbio as liqbio_cmd
from .liqbio import liqbio_prepare as liqbio_prepare_cmd

__author__ = 'dankle'


@click.group()
@click.option('--ref', default='/nfs/PROBIO/autoseq-genome/autoseq-genome.json',
              help='json with reference files to use',
              type=str)
@click.option('--job-params', default=None, help='JSON file specifying various pipeline job ' + \
                                                 'parameters.',
              type=str)
@click.option('--outdir', default='/tmp/autoseq-test', help='output directory', type=click.Path())
@click.option('--libdir', default="/tmp", help="directory to search for libraries")
@click.option('--runner_name', default='shellrunner', help='Runner to use.')
@click.option('--loglevel', default='INFO', help='level of logging')
@click.option('--jobdb', default=None, help="sqlite3 database to write job info and stats")
@click.option('--dot_file', default=None, help="write graph to dot file with this name")
@click.option('--cores', default=1, help="max number of cores to allow jobs to use")
@click.option('--umi', is_flag=True, help="To process the data with UMI- Unique Molecular Identifier")
@click.option('--scratch', default="/tmp", help="scratch dir to use")
@click.pass_context
def cli(ctx, ref, job_params, outdir, libdir, runner_name, loglevel, jobdb, dot_file, cores, umi, scratch):
    setup_logging(loglevel)
    logging.debug("Reading reference data from {}".format(ref))
    ctx.obj = {}
    ctx.obj['refdata'] = load_ref(ref)
    ctx.obj['job_params'] = load_job_params(job_params)
    ctx.obj['outdir'] = outdir
    ctx.obj['libdir'] = libdir
    ctx.obj['pipeline'] = None
    ctx.obj['runner'] = get_runner(runner_name, cores)
    ctx.obj['jobdb'] = jobdb
    ctx.obj['dot_file'] = dot_file
    ctx.obj['cores'] = cores
    ctx.obj['umi'] = umi
    ctx.obj['scratch'] = scratch

    def capture_sigint(sig, frame):
        """
        Capture ctrl-c (or SIGINT sent in other ways).
        1. update remote log
        :param sig:
        :param frame:
        :return:
        """
        try:
            ctx.obj['pipeline'].stop()
            logging.info("Stopping jobs...")
        except AttributeError:
            logging.debug("No pipeline to stop.")

    signal.signal(signal.SIGINT, capture_sigint)
    signal.signal(signal.SIGTERM, capture_sigint)


def load_job_params(job_params_filename):
    # FIXME: Need to validate this and other input files, and clearly define their structure:
    if job_params_filename:
        return json.load(open(job_params_filename))
    else:
        # Return an empty dictionary if no job parameters file is specified, which
        # will result in default job parameters being applied:
        return {}


def convert_to_absolute_path(possible_relative_path, base_path):
    """
    Convert the input potential relative file path to an absolute path by
    prepending the specified base_path, but only if the resulting absolute path points
    to a pre-existing file or directory.

    If the base_path cannot be prepended, then simply return the original input value.

    :param possible_relative_path: A string potentially indicating a relative file/directory path. 
    :param base_path: The base path to prepend.
    :return: Modified path string.
    """

    converted_value = possible_relative_path
    try:
        if not os.path.isabs(possible_relative_path):
            joined_path = os.path.join(base_path, possible_relative_path)
            if os.path.isfile(joined_path) or os.path.isdir(joined_path):
                converted_value = joined_path

    except Exception, e:
        pass

    return converted_value


def make_paths_absolute(input_dict, base_path):
    """Processes the input dictionary, converting relative file paths to absolute
    file paths throughout the dictionary structure.

    Specifically, for each value in the dictionary:
    - If it is also a dictionary, then recursively apply this function,
    replacing the initial dictionary.
    - Otherwise:
    -- If the value is a non-null string that is not already an absolute path,
    then try prepending the specified base_path and see if the resulting file name
    exists, and in that case then replace the string with the resulting absolute path.
    """

    for curr_key, curr_value in input_dict.items():
        if isinstance(curr_value, dict):
            input_dict[curr_key] = make_paths_absolute(curr_value, base_path)
        else:
            converted_value = convert_to_absolute_path(curr_value, base_path)
            input_dict[curr_key] = converted_value

    return input_dict


def load_ref(ref):
    """
    Processes the input genomic reference data JSON file, converting relative file paths
    to absolute paths where required.

    :param ref: Input reference file configuration JSON file.
    :return: Modified reference file dictionary with relative->absolute file path conversions performed.
    """

    basepath = os.path.dirname(ref)
    with open(ref, 'r') as fh:
        refjson = json.load(fh)
        refjson_abs = make_paths_absolute(refjson, basepath)
        return refjson_abs


def get_runner(runner_name, maxcores):
    try:
        module = __import__("pypedream.runners." + runner_name, fromlist="runners")
        runner_class = getattr(module, runner_name.title())

        if runner_name == 'localqrunner':
            runner = runner_class(maxcores)
        else:
            runner = runner_class()
        return runner
    except ImportError:
        print "Couldn't find runner " + runner_name + ". Available Runners:"
        import inspect
        for name, obj in inspect.getmembers(runners):
            if name != "runner" and "runner" in name:
                print "- " + name
        raise ImportError


def setup_logging(loglevel="INFO"):
    """
    Set up logging
    :param loglevel: loglevel to use, one of ERROR, WARNING, DEBUG, INFO (default INFO)
    :return:
    """
    numeric_level = getattr(logging, loglevel.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)
    logging.basicConfig(level=numeric_level,
                        format='%(levelname)s %(asctime)s %(funcName)s - %(message)s')
    logging.info("Started log with loglevel %(loglevel)s" % {"loglevel": loglevel})


cli.add_command(alascca_cmd)
cli.add_command(liqbio_cmd)
cli.add_command(liqbio_prepare_cmd)

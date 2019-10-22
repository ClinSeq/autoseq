"""
General utilities which return static command line strings.
"""

# ## General utilities
from pyjobmanager.job import conditional


def vt_split_and_leftaln(reference_sequence, allow_ref_mismatches=False):
    cmd = "vt decompose -s - " + \
          "|vt normalize " + conditional(allow_ref_mismatches, ' -n ') + \
          " -r {} - ".format(reference_sequence)
    return cmd


def fix_ambiguous_cl(column=4):
    """awk command to replace non-N ambiguous REF bases with N.

    Some callers include these if present in the reference genome but GATK does
    not like them.
    """
    return r"""awk -F$'\t' -v OFS='\t' '{if ($0 !~ /^#/) gsub(/[KMRYSWBVHDX]/, "N", $%s) } {print}'""" % column


def remove_dup_cl():
    """awk command line to remove duplicate alleles where the ref and alt are the same.
    """
    return r""" awk -F$'\t' -v OFS='\t' '$1!~/^#/ && $4 == $5 {next} {print}'"""

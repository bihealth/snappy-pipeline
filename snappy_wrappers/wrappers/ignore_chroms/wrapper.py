"""
The wrapper creates a bed file defining all the genome, except for the ignored chroms (contigs, actually)

For example, if the reference is hs37d5, and the ignored chroms are the decoy, EBV & unplaced/unlocalized contigs,
then the bed file will be::

    1   0   249250621
    2   0   243199373
    3   0   198022430
    ...
    22  0   51304566
    X   0   155270560
    Y   0   59373566
    MT  0   16569


"""

import os
import tempfile
import sys

from snakemake.shell import shell

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.tools.genome_windows import ignore_chroms

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


if not getattr(snakemake, "input", None):
    raise AttributeError("snakemake.input is not defined")
if not getattr(snakemake.input, "reference", None):
    raise AttributeError("Mandatory input 'reference' is not defined")
if not getattr(snakemake, "output", None):
    raise AttributeError("snakemake.output is not defined")
if not getattr(snakemake, "params", None):
    raise AttributeError("snakemake.params is not defined")
    
args = getattr(snakemake.params, "args", {})
with tempfile.NamedTemporaryFile(mode="wt", delete_on_close=False) as f:
    tempfilename = f.name

    for contig, length in ignore_chroms(
        str(snakemake.input.reference),
        set(args.get("ignore_chroms", [])),
        return_ignored=False,
    ):
        print(f"{contig}\t0\t{length}", file=f)

    f.flush()
    f.close()

    cmd = r"""
        bgzip --output {snakemake.output[0]} {tempfilename}
        tabix {snakemake.output[0]}
    """.format(tempfilename=tempfilename, snakemake=snakemake)
    shell(cmd)
import os
import tempfile
import sys

from snakemake.shell import shell

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

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
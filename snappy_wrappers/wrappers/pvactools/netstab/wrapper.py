import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

script = os.path.join(os.path.dirname(__file__), "netMHCstabpan.py")

args = getattr(snakemake.params, "args", {})

lengths = ",".join(map(str, args["lengths"]))
lengths = f"--lengths {lengths}" if lengths else ""

ShellWrapper(snakemake).run(
    r"""
if [[ -s {snakemake.input.epitopes} ]]; then
    export LC_ALL=C.UTF-8

    # netMHCstabpan of the IEDB tools in the pVACtools container; its many small temporary files
    # go to the temporary directory of the job
    python {script} --verbose --workers {snakemake.threads} \
        --tool {args[tool]} {lengths} \
        --base-path $TMPDIR \
        --output {snakemake.output.netstab} \
        {snakemake.input.epitopes}
else
    touch {snakemake.output.netstab}
fi
"""
)

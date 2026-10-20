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
if [[ -s {snakemake.input.epitopes} ]]
then
    # Re-home pVACtools to avoid conflicts with user's setting (.bashrc, ...)
    out=$(mktemp -d)
    out=$(realpath $out)
    mkdir -p $out

    # IMPORTANT: use local /tmp to avoid clobbering cephfs with the creation of many small files
    export TMPDIR=/tmp
    export LC_ALL=C.UTF-8

    python {script} --verbose --workers {snakemake.threads} \
        --tool {args[tool]} {lengths} \
        --base-path $(dirname $out) \
        --container {snakemake.input.container} \
        --output {snakemake.output.netstab} \
        {snakemake.input.epitopes}
else
    touch {snakemake.output.netstab}
fi
"""
)

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

script = os.path.join(os.path.dirname(__file__), "netchop.py")

args = getattr(snakemake.params, "args", {})

# Resource runtime is in minutes, the script wants seconds
timeout = 60 * int(snakemake.resources.runtime)

ShellWrapper(snakemake).run(
    r"""
if [[ -s {snakemake.input.epitopes} ]]
then
    tmpdir=$(dirname {snakemake.output.netchop})/tmp/netchop.$(basename {snakemake.output.netchop})
    mkdir -p $tmpdir

    rm -rf $tmpdir/*
    # IMPORTANT: use local /tmp to avoid clobbering cephfs with the creation of many small files
    export TMPDIR=/tmp
    python {script} --workers {snakemake.threads} \
        --tmpdir $tmpdir --force \
        --method {args[method]} --threshold {args[threshold]} --timeout {timeout} \
        --tool {args[tool]} --netchop {snakemake.input.netchop} \
        --output {snakemake.output.netchop} \
        {snakemake.input.epitopes} \
        {snakemake.input.fasta}
else
    touch {snakemake.output.netchop}
fi
"""
)

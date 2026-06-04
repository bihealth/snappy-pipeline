# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for oxo-G flagging."""

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
mode = args.get("mode", "tag")

ShellWrapper(snakemake).run(
    f"""
set -euo pipefail

set -x

# "Local" TMPDIR as the scripts try to do "rename" across file sests otherwise
export TMPDIR=$(dirname $(dirname {snakemake.output.vcf}))/tmp
mkdir -p $TMPDIR
trap "rm -rf $TMPDIR" EXIT KILL TERM INT HUP

out={snakemake.output.vcf}

dkfzbiasfilter.py \
    --tempFolder $TMPDIR \
    --writeQC \
    {snakemake.input.vcf} \
    {snakemake.input.bam} \
    {snakemake.input.reference} \
    ${{out%.gz}}

if [[ ! -s ${{out%.gz}} ]]; then
    bcftools view --header-only {snakemake.input.vcf} \
    > ${{out%.gz}}
fi

if [[ "{mode}" == "filter" ]]; then
    bcftools filter --exclude 'FILTER ~ "bPcr" || FILTER ~ "bSeq"' -O v -o ${{out%.gz}}.filtered ${{out%.gz}}
    mv ${{out%.gz}}.filtered ${{out%.gz}}
fi

bgzip ${{out%.gz}}
tabix -f {{snakemake.output.vcf}}
"""
)

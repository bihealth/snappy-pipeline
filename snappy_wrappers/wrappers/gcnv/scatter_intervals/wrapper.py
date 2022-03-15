# -*- coding: utf-8 -*-

from snakemake.shell import shell

# Filter interval list file from gCNV model input
interval_list = [
    path
    for path in snakemake.input
    if str(path).endswith(".interval_list")
    and snakemake.wildcards.mapper in str(path)  # Check for correct mapper
    and snakemake.wildcards.library_kit in str(path)  # Check for correct library kit
]

shell(
    r"""
set -x

set -euo pipefail

rm -rf "{snakemake.output}" && mkdir -p "{snakemake.output}"
trap "rm -rf {snakemake.output}" ERR

gatk IntervalListTools \
    --INPUT {interval_list} \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT 5000 \
    --OUTPUT {snakemake.output}
"""
)

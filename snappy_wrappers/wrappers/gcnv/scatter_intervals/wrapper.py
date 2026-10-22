# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args

# Filter interval list file from gCNV model input
interval_list = [path for path in snakemake.input if str(path).endswith(".interval_list")]

ShellWrapper(snakemake).run(
    r"""

set -euo pipefail

rm -rf "{snakemake.output}"
trap "rm -rf {snakemake.output}" ERR

gatk IntervalListTools \
    --INPUT {interval_list} \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT {args[intervals_per_shard]} \
    --OUTPUT {snakemake.output}

# With a single shard (small panels), OUTPUT is the interval list itself, not a directory of
# temp_<i>_of_<n> shards; give it the layout of the shards.
if [[ -f "{snakemake.output}" ]]; then
    mv "{snakemake.output}" "{snakemake.output}.interval_list"
    mkdir -p "{snakemake.output}/temp_0001_of_1"
    mv "{snakemake.output}.interval_list" "{snakemake.output}/temp_0001_of_1/scattered.interval_list"
fi
"""
)

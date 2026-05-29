# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

# Filter interval list file from gCNV model input
interval_list = [path for path in snakemake.input if str(path).endswith(".interval_list")]

ShellWrapper(snakemake).run(
    r"""

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

# -*- coding: utf-8 -*-

import os
from snakemake.shell import shell

shell(
    r"""
set -x

set -euo pipefail

rm -rf "{snakemake.output}" && mkdir -p "{snakemake.output}"
trap "rm -rf {snakemake.output}" ERR

gatk IntervalListTools \
    --INPUT {snakemake.input.interval_list} \
    --SUBDIVISION_MODE INTERVAL_COUNT \
    --SCATTER_CONTENT 5000 \
    --OUTPUT {snakemake.output}
"""
)

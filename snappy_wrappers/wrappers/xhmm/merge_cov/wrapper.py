# -*- coding: utf-8 -*-
# Merge the coverage tracks created in the "gatk_cov" steps.

import os

from snakemake.shell import shell

# Create file with GATK sample_interval_summary output files
# used in `--GATKdepthsList`
work_dir = os.getcwd()
list_path = os.path.join(os.path.dirname(str(snakemake.output)), "gatk_depths.list")
with open(list_path, "w") as f:
    for item in snakemake.input:
        full_path = os.path.join(work_dir, item)
        f.write("%s\n" % full_path)


shell(
    r"""
set -x

xhmm \
    --mergeGATKdepths \
    -o {snakemake.output} \
    --GATKdepthsList {list_path}
"""
)

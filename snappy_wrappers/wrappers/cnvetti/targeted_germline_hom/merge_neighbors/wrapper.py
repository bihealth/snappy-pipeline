# -*- coding: utf-8 -*-
# isort:skip_file
from snappy_wrappers.tools import vcf_cnvetti_coverage_to_hom_del_calls
import os

from snakemake.shell import shell


vcf_cnvetti_coverage_to_hom_del_calls.main([snakemake.output.vcf, snakemake.input.vcf])

shell(
    r"""
set -x
set -euo pipefail

cd $(dirname {snakemake.output.vcf})
cp $(basename {snakemake.output.vcf}) $(basename {snakemake.output.vcf}).bak
tabix -f $(basename {snakemake.output.vcf})

md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}.tbi) >$(basename {snakemake.output.vcf}).tbi.md5
"""
)

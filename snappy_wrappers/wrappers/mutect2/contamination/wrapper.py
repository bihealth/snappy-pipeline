# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

java_options = snakemake.params.args["java_options"]

# Test presence/absence of normal
normal = getattr(snakemake.input, "normal", None)
if normal:
    normal_param = f"--matched-normal {normal}"
else:
    normal_param = ""

ShellWrapper(snakemake).run(
    r"""
# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

gatk --java-options '{java_options}' CalculateContamination \
    --input {snakemake.input.tumor} {normal_param} \
    --tumor-segmentation {snakemake.output.segments} \
    --output {snakemake.output.table} \
    {extra_arguments}
"""
)

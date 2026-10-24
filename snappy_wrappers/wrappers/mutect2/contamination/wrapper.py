# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.input.reference
java_options = snakemake.params.args["java_options"]

# Test presence/absence of normal
normal = getattr(snakemake.input, "normal", None)
if normal:
    normal_param = f"--matched-normal {normal}"
else:
    normal_param = ""

ShellWrapper(snakemake).run(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

out_base=$TMPDIR/$(basename {snakemake.output.table} .contamination.tbl)

gatk --java-options '{java_options}' CalculateContamination \
    --input {snakemake.input.tumor} {normal_param} \
    --tumor-segmentation ${{out_base}}.segments.tbl \
    --output ${{out_base}}.contamination.tbl \
    {extra_arguments}


mv $out_base.* $(dirname {snakemake.output.table})
"""
)

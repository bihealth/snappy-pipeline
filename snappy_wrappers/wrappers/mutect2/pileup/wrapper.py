# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GetPileupSummaries: Snakemake wrapper.py"""

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference
common_variants = snakemake.input.common_variants

if java_options := args.get("java_options", None):
    java_options = f"--java-options '{java_options}'"

extra_arguments = " ".join(args.get("extra_arguments", []))

shell.executable("/bin/bash")

ShellWrapper(snakemake).run(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

out_base=$TMPDIR/$(basename {snakemake.output.pileup} .pileup)

gatk {java_options} GetPileupSummaries \
    --input {snakemake.input.bam} \
    --reference {reference} \
    --variant {common_variants} \
    --intervals {common_variants} \
    --output $out_base.pileup \
    {extra_arguments}

pushd $TMPDIR && \
    for f in $out_base.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv $out_base.* $(dirname {snakemake.output.pileup})
"""
)


# -*- coding: utf-8 -*-
"""Wrapper for running bcftools call (with samtools mpileup) in joint calling mode.

This can be used for performing variant calling on all samples from a given donor for somatic
variant calling.
"""

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
args_ignore_chroms = ""
if ignore_chroms := args.get("ignore_chroms"):
    args_ignore_chroms = " ".join(["--ignore-chroms"] + ignore_chroms)

reference_path = snakemake.input.reference
max_depth = args["max_depth"]
max_indel_depth = args["max_indel_depth"]
window_length = args["window_length"]
num_threads = args["num_threads"]

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------


export REF={reference_path}

bcftools_joint()
{{
    samtools mpileup \
        --BCF \
        --fasta-ref $REF \
        --max-depth {max_depth} \
        --max-idepth {max_indel_depth} \
        --output-tags DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
        --per-sample-mF \
        --redo-BAQ \
        --uncompressed \
        -r $1 \
        {snakemake.input.bam} \
    | bcftools call \
        --format-fields GQ,GP \
        -m \
        -v
}}
export -f bcftools_joint

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

snappy-genome_windows \
    --fai-file $REF.fai \
    --window-size {window_length} \
    {args_ignore_chroms} \
| parallel \
    --plain \
    --keep-order \
    --verbose \
    --max-procs {num_threads} \
    bcftools_joint \
| snappy-vcf_first_header \
| bcftools norm \
    --fasta-ref $REF \
    --multiallelics -any \
| snappy-vcf_sort $REF.fai \
| bgzip -c \
> {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

# -*- coding: utf-8 -*-
"""Wrapper for running bcftools call (with samtools mpileup) in joint calling mode.

This can be used for performing variant calling on all samples from a given donor for somatic
variant calling.
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args_ignore_chroms = ""
if snakemake.params.args["ignore_chroms"]:
    args_ignore_chroms = " ".join(["--ignore-chroms"] + snakemake.params.args["ignore_chroms"])

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------


export REF={snakemake.config[static_data_config][reference][path]}

bcftools_joint()
{{
    samtools mpileup \
        --BCF \
        --fasta-ref $REF \
        --max-depth {snakemake.config[step_config][somatic_variant_calling][bcftools_joint][max_depth]} \
        --max-idepth {snakemake.config[step_config][somatic_variant_calling][bcftools_joint][max_indel_depth]} \
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
    --window-size {snakemake.config[step_config][somatic_variant_calling][bcftools_joint][window_length]} \
    {args_ignore_chroms} \
| parallel \
    --plain \
    --keep-order \
    --verbose \
    --max-procs {snakemake.config[step_config][somatic_variant_calling][bcftools_joint][num_threads]} \
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

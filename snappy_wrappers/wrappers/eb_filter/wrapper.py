# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for EBFilter flagging.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

if snakemake.params["args"]["interval"]:
    cmd_fetch = "tabix --print-header {} {}".format(
        snakemake.input.vcf, snakemake.params["args"]["interval"]
    )
else:
    cmd_fetch = "zcat {}".format(snakemake.input.vcf)

shell(
    r"""
set -x

export TMPDIR=$HOME/scratch/tmp
mkdir -p $TMPDIR

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

export REF={snakemake.config[static_data_config][reference][path]}

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Used to be:
# filter='FILTER == "germline_risk" || FILTER == "t_lod_fstar" || FILTER == "OffExome" || ANN ~ "stream_gene_variant"'
if [[ {snakemake.input.vcf} == *"mutect2"* ]]; then
    filter='FILTER == "germline" || FILTER == "weak_evidence" || FILTER == "OffExome" || ANN ~ "stream_gene_variant"'
    {cmd_fetch} \
    | bcftools view \
        -e "$filter" \
        -O z \
        -o $TMPDIR/for_eb_filter.vcf.gz

    {cmd_fetch} \
    | bcftools view \
        -i "$filter" \
        -O z \
        -o $TMPDIR/not_for_eb_filter.vcf.gz
else
    {cmd_fetch} | bgzip -c > $TMPDIR/for_eb_filter.vcf.gz
    zgrep "^#" {snakemake.input.vcf} | bgzip -c > $TMPDIR/not_for_eb_filter.vcf.gz
fi

set +e
lines=$(zgrep -v '^#' $TMPDIR/for_eb_filter.vcf.gz | wc -l)
set -e

# EBFilter does not like empty files ...
if [[ $lines -gt 0 ]]; then
    EBFilter \
        -f vcf \
        -t 1 \
        -q {snakemake.config[step_config][somatic_variant_filtration][eb_filter][min_mapq]} \
        -Q {snakemake.config[step_config][somatic_variant_filtration][eb_filter][min_baseq]} \
        $TMPDIR/for_eb_filter.vcf.gz \
        {snakemake.input.bam} \
        {snakemake.input.txt} \
        $TMPDIR/after_eb_filter.vcf
else
    zcat $TMPDIR/for_eb_filter.vcf.gz \
    > $TMPDIR/after_eb_filter.vcf
fi

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

bcftools concat \
    $TMPDIR/after_eb_filter.vcf \
    $TMPDIR/not_for_eb_filter.vcf.gz \
| snappy-vcf_sort $REF.fai \
| bgzip -c \
> {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5 && \
    popd
"""
)

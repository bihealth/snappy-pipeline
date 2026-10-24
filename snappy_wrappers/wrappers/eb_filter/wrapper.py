# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for EBFilter flagging."""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
filter_name = args["filter_name"] if "filter_name" in args else ""
mode = args.get("mode", "tag")
has_annotation = str(args["has_annotation"] if "has_annotation" in args else False)
skip_filters = " ".join(args["skip_filters"])
skip_annotation = args["skip_annotation"]

if "interval" in args:
    cmd_fetch = "tabix --print-header {} {}".format(snakemake.input.vcf, args["interval"])
else:
    cmd_fetch = "zcat {}".format(snakemake.input.vcf)

ShellWrapper(snakemake).run(
    r"""
set -x

export TMPDIR=$HOME/scratch/tmp
mkdir -p $TMPDIR

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

export REF={snakemake.input.reference}

# Used to be:
# filter='FILTER == "germline_risk" || FILTER == "t_lod_fstar" || FILTER == "OffExome" || ANN ~ "stream_gene_variant"'
# Mutect2 calls are recognized by the FILTER IDs in their header (the paths do not name the caller).
# Variants with these filters, and those in or near genes with annotation, skip EBFilter.
filter=""
for f in $(echo "{skip_filters}" | tr ' ' '\n')
do
    use="Yes"
    zgrep -q "^##FILTER=<ID=$f," {snakemake.input.vcf} || use="No"
    if [[ $use = "Yes" ]]
    then
        filter="$filter || FILTER == \"$f\""
    fi
done
filter=$(echo "$filter" | sed -e "s/^ || //")
if [[ -n "$filter" ]]; then
    ann="CSQ"
    zgrep -q "^##INFO=<ID=CSQ," {snakemake.input.vcf} || ann="ANN"
    if [[ '{has_annotation}' == "True" && -n "{skip_annotation}" ]]; then
        filter="$filter || $ann ~ \"{skip_annotation}\""
    else
        filter="$filter"
    fi
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
        -t {snakemake.threads} \
        -q {args[min_mapq]} \
        -Q {args[min_baseq]} \
        {extra_arguments} \
        $TMPDIR/for_eb_filter.vcf.gz \
        {snakemake.input.bam} \
        {snakemake.input.txt} \
        $TMPDIR/after_running_eb_filter.vcf
    if [[ "{mode}" == "tag" ]]
    then
        if [[ -n "{filter_name}" ]]
        then
            bcftools filter --soft-filter {filter_name} --mode + \
                --exclude "INFO/EB < {args[ebfilter_threshold]}" \
                -O z -o $TMPDIR/after_eb_filter.vcf.gz \
                $TMPDIR/after_running_eb_filter.vcf
        else
            mv $TMPDIR/after_running_eb_filter.vcf $TMPDIR/after_eb_filter.vcf
            bgzip $TMPDIR/after_eb_filter.vcf
        fi
    else
        bcftools filter \
            --exclude "INFO/EB < {args[ebfilter_threshold]}" \
            -O z -o $TMPDIR/after_eb_filter.vcf.gz \
            $TMPDIR/after_running_eb_filter.vcf
    fi
else
    mv $TMPDIR/for_eb_filter.vcf.gz $TMPDIR/after_eb_filter.vcf.gz
fi

bcftools index --force $TMPDIR/after_eb_filter.vcf.gz
bcftools index --force $TMPDIR/not_for_eb_filter.vcf.gz

bcftools concat \
    --allow-overlaps \
    $TMPDIR/after_eb_filter.vcf.gz \
    $TMPDIR/not_for_eb_filter.vcf.gz \
| bcftools sort --output {snakemake.output.vcf} --output-type z

tabix -f {snakemake.output.vcf}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for EBFilter flagging.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

if "args" in snakemake.params and snakemake.params["args"].get("interval", None):
    cmd_fetch = "tabix --print-header {} {}".format(
        snakemake.input.vcf, snakemake.params["args"]["interval"]
    )
else:
    cmd_fetch = "zcat {}".format(snakemake.input.vcf)

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]
has_annotation = str(config.get("has_annotation", False))

if "ebfilter_threshold" in config:
    threshold = config.get("ebfilter_threshold", 0)
elif "eb_filter" in config and "ebfilter_threshold" in config["eb_filter"]:
    threshold = config["eb_filter"].get("ebfilter_threshold", 0)
elif "ebfilter" in config and "ebfilter_threshold" in config["ebfilter"]:
    threshold = config["ebfilter"].get("ebfilter_threshold", 0)
else:
    try:
        filter_nb = int(snakemake.wildcards["filter_nb"]) - 1
        threshold = config["filter_list"][filter_nb]["ebfilter"].get("ebfilter_threshold", 0)
    except:
        threshold = 0

if "filter_name" in config:
    filter_name = config.get("filter_name", "")
elif "eb_filter" in config and "filter_name" in config["eb_filter"]:
    filter_name = config["eb_filter"].get("filter_name", "")
elif "ebfilter" in config and "filter_name" in config["ebfilter"]:
    filter_name = config["ebfilter"].get("filter_name", "")
else:
    try:
        filter_name = "ebfilter_{}".format(int(snakemake.wildcards["filter_nb"]))
    except:
        filter_name = "+"

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
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

# Used to be:
# filter='FILTER == "germline_risk" || FILTER == "t_lod_fstar" || FILTER == "OffExome" || ANN ~ "stream_gene_variant"'
if [[ {snakemake.input.vcf} == *"mutect2"* ]]; then
    filter=""
    for f in $(echo "germline weak_evidence OffExome" | tr ' ' '\n')
    do
        use="Yes"
        zgrep -q "^##FILTER=<ID=$f," {snakemake.input.vcf} || use="No"
        if [[ $use = "Yes" ]]
        then
            filter="$filter || FILTER == \"$f\""
        fi
    done
    filter=$(echo "$filter" | sed -e "s/^ || //")
    ann="CSQ"
    zgrep -q "^##INFO=<ID=CSQ," {snakemake.input.vcf} || ann="ANN"
    if [[ '{has_annotation}' == "True" ]]; then
        filter="$filter || $ann ~ \"stream_gene_variant\""
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
        -t 1 \
        -q {snakemake.config[step_config][somatic_variant_filtration][eb_filter][min_mapq]} \
        -Q {snakemake.config[step_config][somatic_variant_filtration][eb_filter][min_baseq]} \
        $TMPDIR/for_eb_filter.vcf.gz \
        {snakemake.input.bam} \
        {snakemake.input.txt} \
        $TMPDIR/after_running_eb_filter.vcf
    bcftools filter --soft-filter {filter_name} --mode + \
        --exclude "INFO/EB < {threshold}" \
        -O z -o $TMPDIR/after_eb_filter.vcf.gz \
        $TMPDIR/after_running_eb_filter.vcf
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

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5 && \
    popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

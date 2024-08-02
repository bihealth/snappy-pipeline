# -*- coding: utf-8 -*-
"""Wrapper for running Varscan calling."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

# The step and caller key to use.
step_key = snakemake.params.step_key
caller_key = snakemake.params.caller_key

args_ignore_chroms = ""
if snakemake.params.args["ignore_chroms"]:
    args_ignore_chroms = " ".join(["--ignore-chroms"] + snakemake.params.args["ignore_chroms"])

if snakemake.config["step_config"][step_key][caller_key]["no_baq"]:
    arg_no_baq = "-B"
else:
    arg_no_baq = ""

min_bq = snakemake.config["step_config"][step_key][caller_key]["min_bq"]
arg_min_bq = "--min-BQ {}".format(min_bq)
arg_min_avg_qual = "--min-avg-qual {}".format(
    snakemake.config["step_config"][step_key][caller_key]["min_bq"]
)

# If intervals are emtpy then will be generated for whole genome.
intervals = " ".join(snakemake.params["args"].get("intervals", []))

sample_list = " ".join(snakemake.params["args"]["sample_list"])

shell(
    r"""
export TMPDIR=$HOME/scratch/tmp

# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

export REF={snakemake.config[static_data_config][reference][path]}

# Java fun!
export MALLOC_ARENA_MAX=4

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

echo {snakemake.params[args][sample_list]} \
| tr ' ' '\n' \
> $TMPDIR/samples.txt

varscan_joint()
{{
    set -x
    vartype=$1
    region=$2

    # We need to disable pipefail here because varscan fails if samtools mpileup generates
    # an empty output file. We rescue this case with snappy-fix_vcf.
    set +o pipefail

    # The pipeline below is:
    #
    # - samtools mpileup
    # - varscan mpileup2<snp or indel>
    # - replace non-ACGTN characters in reference with "N"
    samtools mpileup \
        {arg_no_baq} \
        {arg_min_bq} \
        -r $region \
        --reference $REF \
        --max-depth {snakemake.config[step_config][somatic_variant_calling][varscan_joint][max_depth]} \
        --max-idepth {snakemake.config[step_config][somatic_variant_calling][varscan_joint][max_indel_depth]} \
        $(echo {snakemake.input} | tr ' ' '\n' | grep 'bam$') \
    | varscan mpileup2$vartype \
        --output-vcf \
        {arg_min_avg_qual} \
        --min-coverage {snakemake.config[step_config][somatic_variant_calling][varscan_joint][min_coverage]} \
        --min-reads2 {snakemake.config[step_config][somatic_variant_calling][varscan_joint][min_reads2]} \
        --min-avg-qual {snakemake.config[step_config][somatic_variant_calling][varscan_joint][min_avg_qual]} \
        --min-var-freq {snakemake.config[step_config][somatic_variant_calling][varscan_joint][min_var_freq]} \
        --min-freq-for-hom {snakemake.config[step_config][somatic_variant_calling][varscan_joint][min_freq_for_hom]} \
        --p-value {snakemake.config[step_config][somatic_variant_calling][varscan_joint][p_value]} \
        --vcf-sample-list $TMPDIR/samples.txt \
        /dev/stdin \
    | awk -F $'\t' 'BEGIN {{ OFS=FS; }} /^#/ {{ print; }} !/^#/ {{ gsub(/[^ACGTN]/, "N", $4); print; }}' \
    | snappy-fix_vcf --faidx $REF.fai --input-vcf /dev/stdin --sample {sample_list}
}}
export -f varscan_joint

generate_intervals()
{{
    set -x

    if [[ -n "{intervals}" ]]; then
        echo {intervals} \
        | tr ' ' '\n'
    else
        snappy-genome_windows \
            --fai-file $REF.fai \
            --window-size {snakemake.config[step_config][somatic_variant_calling][varscan_joint][window_length]} \
            {args_ignore_chroms}
    fi
}}


for vartype in snp indel; do
    generate_intervals \
    | parallel \
        --plain \
        --keep-order \
        --verbose \
        --max-procs {snakemake.config[step_config][somatic_variant_calling][varscan_joint][num_cores]} \
        "varscan_joint $vartype {{}}" \
    | snappy-vcf_first_header \
    | bcftools norm \
        --fasta-ref $REF \
        -d both \
        /dev/stdin \
    | snappy-vcf_sort $REF.fai \
    | bgzip -c \
    > $TMPDIR/$vartype.vcf.gz
    tabix -f $TMPDIR/$vartype.vcf.gz
done

bcftools concat \
    --allow-overlaps \
    -O u \
    $TMPDIR/snp.vcf.gz \
    $TMPDIR/indel.vcf.gz \
| bcftools norm \
    -O z \
    -o {snakemake.output.vcf} \
    --fasta-ref $REF \
    -d both \
    /dev/stdin

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

DEF_HELPER_FUNCS = r"""
compute-md5()
{
    if [[ $# -ne 2 ]]; then
        >&2 echo "Invalid number of arguments: $#"
        exit 1
    fi
    md5sum $1 \
    | awk '{ gsub(/.*\//, "", $2); print; }' \
    > $2
}
"""
args = snakemake.params["args"]

reference_path = snakemake.input.reference
reference_index_path = snakemake.input.reference_index

assembly = args["assembly"]

ignore_chroms = args["ignore_chroms"]

max_depth = args["max_depth"]
max_indel_depth = args["max_indel_depth"],

gatk4_hc_joint_window_length = args["gatk4_hc_joint_window_length"]
gatk4_hc_joint_num_threads = args["gatk4_hc_joint_num_threads"]

shell(
    r"""
set -x

# Write files for reproducibility -----------------------------------------------------------------

{DEF_HELPER_FUNCS}

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
compute-md5 {snakemake.log.conda_list} {snakemake.log.conda_list_md5}
compute-md5 {snakemake.log.conda_info} {snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
compute-md5 {snakemake.log.wrapper} {snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
compute-md5 {snakemake.log.env_yaml} {snakemake.log.env_yaml_md5}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Run actual tools --------------------------------------------------------------------------------

# Create binning of the reference into windows of roughly the same size.
gatk PreprocessIntervals \
    --reference {ref_path} \
    --bin-length {gatk4_hc_joint_window_length} \
    --output $TMPDIR/raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for ignore_chrom in {ignore_chroms}; do \
        awk "(\$1 ~ /$ignore_chrom/) {{ printf(\"--exclude-intervals %s:1-%d\\n\", \$1, \$2) }}" \
            {reference_index_path}; \
    done)

# Postprocess the Picard-style interval list into properly padded interval strings suitable for
# passing to ``--intervals``.
awk -v PADDING=1000 '
    (!/^@/) {{
        chrom=$1;
        start=$2;
        if (start > PADDING) {{
            start = start - PADDING
        }} else {{
            start = 1
        }}
        end=$3;
        printf("%s:%d-%d\n", chrom, start, end);
    }}
' $TMPDIR/raw.interval_list \
> $TMPDIR/final_intervals.txt
wc -l $TMPDIR/final_intervals.txt

# Create per-shard output directory
mkdir -p $TMPDIR/shards-output

# Function to run "bcftools call"
run-shard()
{{
    set -x

    job_no=$1
    interval=$2

    bcftools mpileup \
        -Ou \
        --annotate FORMAT/AD,FORMAT/DP \
        -f {reference_path} \
        --per-sample-mF \
        --max-depth {max_depth} \
        --max-idepth {max_indel_depth} \
        --redo-BAQ \
        --regions $2 \
        {snakemake.input.bam} \
    | bcftools call \
        --multiallelic-caller \
        $(if [[ "{assembly}" != "unknown" ]]; then \
            echo --ploidy={assembly}; \
            echo --samples-file {snakemake.input.ped}; \
        fi) \
        --variants-only \
        -Oz \
        -o $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz
    tabix -f $TMPDIR/shards-output/$(printf %06d $job_no).vcf.gz
}}
export -f run-shard

# Perform parallel execution
(set -x; sleep $(echo "scale=3; $RANDOM/32767*10" | bc)s) # sleep up to 10s to work around bug
num_threads={gatk4_hc_joint_num_threads}
cat $TMPDIR/final_intervals.txt \
| parallel --plain -j $num_threads 'run-shard {{#}} {{}}'

# Merge the individual shards' output VCF
bcftools concat \
    --allow-overlaps  \
    -d none \
    -O u \
    $TMPDIR/shards-output/*.vcf.gz \
| bcftools sort \
    -T $TMPDIR/bcftools.sort.XXXXXX \
    -O u \
    /dev/stdin \
| bcftools norm \
    -d exact \
    -f {reference_path} \
    -O z \
    -o {snakemake.output.vcf}
tabix {snakemake.output.vcf}

# Compute MD5 sums on output files
compute-md5 {snakemake.output.vcf} {snakemake.output.vcf_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
compute-md5 {snakemake.output.vcf_tbi} {snakemake.output.vcf_tbi_md5}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
{DEF_HELPER_FUNCS}

sleep 1s  # try to wait for log file flush
compute-md5 {snakemake.log.log} {snakemake.log.log_md5}
"""
)

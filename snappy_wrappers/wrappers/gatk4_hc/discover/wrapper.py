from snakemake import shell

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
    --reference {snakemake.config[static_data_config][reference][path]} \
    --bin-length {snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][window_length]} \
    --output $TMPDIR/raw.interval_list \
    --interval-merging-rule OVERLAPPING_ONLY \
    $(for ignore_chrom in {snakemake.config[step_config][variant_calling][ignore_chroms]}; do \
        awk "(\$1 ~ /$ignore_chrom/) {{ printf(\"--exclude-intervals %s:1-%d\\n\", \$1, \$2) }}" \
            {snakemake.config[static_data_config][reference][path]}.fai; \
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

# Function to run HaplotypeCaller in to generate gVCF
run-shard()
{{
    job_no=$1
    interval=$2

    # Create per-sample gVCF file
    GATK_JAVA_MEMORY=3750m
    gatk \
        HaplotypeCaller \
        --java-options "-Xmx$GATK_JAVA_MEMORY -Djava.io.tmpdir=$TMPDIR" \
        --tmp-dir $TMPDIR \
        --output $TMPDIR/shards-output/$(printf %06d $job_no).g.vcf.gz \
        --reference {snakemake.config[static_data_config][reference][path]} \
        --dbsnp {snakemake.config[static_data_config][dbsnp][path]} \
        --intervals $interval \
        --input {snakemake.input.bam} \
        -G StandardAnnotation \
        -G AS_StandardAnnotation \
        -G StandardHCAnnotation \
        -G AlleleSpecificAnnotation \
        $(if [[ {snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][allow_seq_dict_incompatibility]} == "True" ]]; then \
            echo --disable-sequence-dictionary-validation true; \
        fi) \
        -ERC GVCF
}}
export -f run-shard

# Perform parallel execution
(set -x; sleep $(echo "scale=3; $RANDOM/32767*10" | bc)s) # sleep up to 10s to work around bug
num_threads={snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][num_threads]}
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
    -f {snakemake.config[static_data_config][reference][path]} \
    -O z \
    -o {snakemake.output.gvcf}
tabix {snakemake.output.gvcf}

# Compute MD5 sums on output files
compute-md5 {snakemake.output.gvcf} {snakemake.output.gvcf_md5}
compute-md5 {snakemake.output.gvcf_tbi} {snakemake.output.gvcf_tbi_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
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

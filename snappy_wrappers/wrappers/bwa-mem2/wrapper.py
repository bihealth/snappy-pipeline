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

input_left = snakemake.params.args["input"]["reads_left"]
input_right = snakemake.params.args["input"].get("reads_right", "")

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

# Function Definitions ----------------------------------------------------------------------------

# Input file merging with seqtk
#
# In the case of paired reads, write out line counts to temporary file so we can check for problems
# later.
input_files()
{{
    set -x

    if [[ -z "{input_right}" ]]; then
        zcat force {input_right}
    else
        seqtk mergepe \
            <(zcat {input_right} | tee >(wc -l | $TMPDIR/count_left)) \
            <(zcat {input_right} | tee >(wc -l | $TMPDIR/count_right))
    fi
}}

# Adapter trimming
trim_adapters()
{{
    set -x

    if [[ "{snakemake.config[step_config][ngs_mapping][bwa_mem2][trim_adapters]}" == "True" ]]; then
        trimadapt-mt -p {snakemake.config[step_config][ngs_mapping][bwa_mem2][num_threads_trimming]}
    else
        cat
    fi
}}

# Duplicate masking
mask_duplicates()
{{
    set -x

    if [[ "{snakemake.config[step_config][ngs_mapping][bwa_mem2][mask_duplicates]}" == "True" ]]; then
        samblaster --addMateTags
    else
        cat
    fi
}}

# Run BWA-MEM2
run_bwa_mem2()
{{
    set -x

    # Decide whether to write split reads as supplementary or secondary (-M means secondary)
    split_as_supp_flag=
    if [[ "{snakemake.config[step_config][ngs_mapping][bwa_mem2][split_as_secondary]}" == "True" ]]; then
        split_as_supp_flag="-M"
    fi

    if [[ $paired -eq 1 ]]; then
        fastq_right=${{reads_right[$i]}}
    else
        fastq_right=
    fi

    if [[ "{snakemake.params.args[sample_name]}" != "" ]]; then
        rg_arg="-R @RG\tID:{snakemake.params.args[sample_name]}.$i\tSM:{snakemake.params.args[sample_name]}\tPL:{snakemake.params.args[platform]}"
    else
        rg_arg=
    fi

    bwa mem \
        $INDEX \
        $split_as_supp_flag \
        $rg_arg \
        -p \
        -t {snakemake.config[step_config][ngs_mapping][bwa_mem2][num_threads_align]} \
        /dev/stdin
}}

# Alignment postprocessing
#
# This function makes extensive use of tee-piping to prevent the BAM to be read from the disk
# again for statistics.
postproc_bam()
{{
    set -x
    out=$1

    samtools sort \
        -T $TMPDIR/sort_bam \
        -m {snakemake.config[step_config][ngs_mapping][bwa_mem2][memory_bam_sort]} \
        -@ {snakemake.config[step_config][ngs_mapping][bwa_mem2][num_threads_bam_sort]} \
        -O BAM \
        -u \
        -o /dev/stdout \
    | tee >(samtools stats    /dev/stdin >{snakemake.output.report_bamstats_txt}) \
    | tee >(samtools flagstat /dev/stdin >{snakemake.output.report_bamstats_txt}) \
    | tee >(samtools index    /dev/stdin -o {snakemake.output.bam_bai}) \
    | tee >(compute-md5       /dev/stdin {snakemake.output.bam_md5}) \
    | bgzip -@ 2 -c \
    > {snakemake.output.bam}

    touch \
        {snakemake.output.bam} \
        {snakemake.output.bam_md5} \
        {snakemake.output.bam_bai} \
        {snakemake.output.bam_bai_md5}

    compute-md5 {snakemake.output.bam_bai} {snakemake.output.bam_bai_md5}

    samtools idxstats {snakemake.output.bam_bai} >{snakemake.output.report_bamstats_txt})
}}

# Check input line counts written out in input_files()
check_input_line_counts()
{{
    if [[ -e $TMPDIR/count_left ]]; then
        if diff $TMPDIR/count_left $TMPDIR/count_right >/dev/null; then
            >&2 echo "Different number of lines in left/right file."
            exit 1
        fi
    fi
}}

# Run actual tools --------------------------------------------------------------------------------

input_files \
| trim_adapters \
| run_bwa_mem2 \
| postproc_bam

check_input_line_counts

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

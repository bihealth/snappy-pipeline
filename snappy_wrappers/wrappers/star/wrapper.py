from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

out_gc = snakemake.output.get("gene_counts", "__dummy__")
out_sj = snakemake.output.get("junctions", "__dummy__")
out_tx = snakemake.output.get("transcriptome", "__dummy__")

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

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = snakemake.params.args["input"]["reads_left"]
reads_right = snakemake.params.args["input"].get("reads_right", "")

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

# Some More Preparation ---------------------------------------------------------------------------

mkdir -p $TMPDIR/tmp.d $TMPDIR/pre.d

# Define some global shortcuts
INDEX={snakemake.config[step_config][ngs_mapping][star][path_index]}

# Define left and right reads as Bash arrays
declare -a reads_left=({reads_left})
# declared but never used
declare -a reads_right=({reads_right})

# Function Definitions ----------------------------------------------------------------------------

# Duplicate masking from SAM file in STDIN to STDOUT
mask_duplicates()
{{
    set -x

    samblaster --addMateTags 
}}

# Sort by coordinate input SAM from STDIN
sort_by_coord()
{{
    set -x

    samtools sort -O BAM \
        -m {snakemake.config[step_config][ngs_mapping][star][memory_bam_sort]} \
        -@ {snakemake.config[step_config][ngs_mapping][star][num_threads_bam_sort]} \
        -
}}

# Merge and index STAR output
index_bam()
{{
    set -x

    out_bam=$1

    # Index resulting BAM file
    samtools index $out_bam

    # Build MD5 files
    compute-md5 $out_bam $out_bam.md5
    compute-md5 $out_bam.bai $out_bam.bai.md5
}}

# Function for running STAR
run_star()
{{
    set -x

    rg_args=""

    if [[ "{snakemake.params.args[sample_name]}" != "" ]]; then
        rg_args="--outSAMattrRGline "
        for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
            if [[ $i -gt 0 ]]; then
                rg_args="${{rg_args}} , "
            fi
            rg_arg="\"ID:{snakemake.params.args[sample_name]}.$i\" \"SM:{snakemake.params.args[sample_name]}\" \"PL:{snakemake.params.args[platform]}\""
            rg_args="${{rg_args}}${{rg_arg}}"
        done
    fi

    left_files=$(IFS="," ; echo "${{reads_left[*]}}")

    right_files=""
    if [[ "{reads_right}" != "" ]]; then
        right_files=$(IFS="," ; echo "${{reads_right[*]}}")
    fi

    trim_cmd=""
    if [[ "{snakemake.config[step_config][ngs_mapping][star][trim_adapters]}" == "True" ]]; then
        trim_cmd="\"trimadap-mt -p {snakemake.config[step_config][ngs_mapping][star][num_threads_trimming]}\""
    else
        trim_cmd="zcat"
    fi

    quant_mode=""
    if [[ -n "{snakemake.config[static_data_config][features][path]}" ]]
    then
        quant_mode="$quant_mode GeneCounts"
    fi
    if [[ "{snakemake.config[step_config][ngs_mapping][star][transcriptome]}" = "True" ]]
    then
        quant_mode="$quant_mode TranscriptomeSAM"
    fi

    STAR \
        --readFilesIn ${{left_files}} ${{right_files}} \
        {snakemake.config[step_config][ngs_mapping][star][raw_star_options]} \
        $rg_args \
        --readFilesCommand ${{trim_cmd}} \
        --alignIntronMax {snakemake.config[step_config][ngs_mapping][star][align_intron_max]} \
        --alignIntronMin {snakemake.config[step_config][ngs_mapping][star][align_intron_min]} \
        --alignMatesGapMax {snakemake.config[step_config][ngs_mapping][star][align_mates_gap_max]} \
        --alignSJDBoverhangMin {snakemake.config[step_config][ngs_mapping][star][align_sjdb_overhang_min]} \
        --alignSJoverhangMin {snakemake.config[step_config][ngs_mapping][star][align_sj_overhang_min]} \
        --genomeDir {snakemake.config[step_config][ngs_mapping][star][path_index]} \
        --genomeLoad {snakemake.config[step_config][ngs_mapping][star][genome_load]} \
        --outFileNamePrefix $TMPDIR/pre.d/out. \
        --outFilterIntronMotifs {snakemake.config[step_config][ngs_mapping][star][out_filter_intron_motifs]} \
        --outFilterMismatchNmax {snakemake.config[step_config][ngs_mapping][star][out_filter_mismatch_n_max]} \
        --outFilterMismatchNoverLmax {snakemake.config[step_config][ngs_mapping][star][out_filter_mismatch_n_over_l_max]} \
        --outFilterMultimapNmax {snakemake.config[step_config][ngs_mapping][star][out_filter_multimap_n_max]} \
        --outFilterType {snakemake.config[step_config][ngs_mapping][star][out_filter_type]} \
        --outSAMstrandField {snakemake.config[step_config][ngs_mapping][star][out_sam_strand_field]} \
        --outSAMunmapped $(if [[ "{snakemake.config[step_config][ngs_mapping][star][include_unmapped]}" == "True" ]]; then \
                echo "Within"; \
            else
                echo "None"; \
            fi) \
        $(if [[ -n "$quant_mode" ]]; then \
            echo "--quantMode $quant_mode"
        fi) \
        $(if [[ -n "{snakemake.config[static_data_config][features][path]}" ]]; then \
            echo --sjdbGTFfile "{snakemake.config[static_data_config][features][path]}"
        fi) \
        $(if [[ "{snakemake.config[step_config][ngs_mapping][star][mask_duplicates]}" == "True" ]]; then \
            echo " --outStd SAM " ; \
        else
            echo " --outSAMtype BAM SortedByCoordinate "; \
        fi) \
        --runThreadN {snakemake.config[step_config][ngs_mapping][star][num_threads_align]} 

    >&2 ls -lhR $TMPDIR
}}

# Perform Alignment -------------------------------------------------------------------------------

# Run STAR
if [[ "{snakemake.config[step_config][ngs_mapping][star][mask_duplicates]}" == "True" ]]; then
    run_star | mask_duplicates | sort_by_coord > {snakemake.output.bam}
else
    run_star
    mv $TMPDIR/pre.d/out.Aligned.sortedByCoord.out.bam {snakemake.output.bam}
fi

index_bam {snakemake.output.bam}

mv $TMPDIR/pre.d/out.ReadsPerGene.out.tab {out_gc}
compute-md5 {out_gc} {out_gc}.md5
mv $TMPDIR/pre.d/out.SJ.out.tab {out_sj}
compute-md5 {out_sj} {out_sj}.md5

# Optional output: mapping on transcriptome -------------------------------------------------------

if [[ "{snakemake.config[step_config][ngs_mapping][star][transcriptome]}" = "True" ]]; then
    if [[ "{snakemake.config[step_config][ngs_mapping][star][mask_duplicates]}" == "True" ]]; then
        samtools view -h -S $TMPDIR/pre.d/out.Aligned.toTranscriptome.out.bam | mask_duplicates | samtools view -h -b - > {out_tx}
    else
        mv $TMPDIR/pre.d/out.Aligned.toTranscriptome.out.bam {out_tx}
    fi
    compute-md5 {out_tx} {out_tx}.md5
fi

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

# Build MD5 files for the reports
compute-md5 {snakemake.output.report_bamstats_txt}  {snakemake.output.report_bamstats_txt_md5}
compute-md5 {snakemake.output.report_flagstats_txt} {snakemake.output.report_flagstats_txt_md5}
compute-md5 {snakemake.output.report_idxstats_txt}  {snakemake.output.report_idxstats_txt_md5}

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

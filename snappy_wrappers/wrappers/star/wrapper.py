# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for STAR: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = snakemake.params.args["input"]["reads_left"]
reads_right = snakemake.params.args["input"].get("reads_right", "")

shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

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

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d $TMPDIR/pre.d

# Define some global shortcuts
INDEX={snakemake.config[step_config][ngs_mapping][star][path_index]}

# Define left and right reads as Bash arrays
declare -a reads_left=({reads_left})
# declared but never used
declare -a reads_right=({reads_right})

# Function Definitions ----------------------------------------------------------------------------

# Duplicate masking
mask_duplicates()
{{
    set -x

    if [[ "{snakemake.config[step_config][ngs_mapping][star][mask_duplicates]}" == "True" ]]; then
        samtools view \
            -Sh \
            -@ {snakemake.config[step_config][ngs_mapping][star][num_threads_bam_view]} \
            $1 \
        | samblaster --addMateTags \
        | samtools view \
            -u \
            -Sb \
            -@ {snakemake.config[step_config][ngs_mapping][star][num_threads_bam_view]}
    else
        cat $1 # TODO: can we somehow remove this?
    fi
}}

# Alignment postprocessing
postproc_bam()
{{
    set -x
    in=$1
    out=$2

    mask_duplicates $in \
    | samtools sort \
        -m {snakemake.config[step_config][ngs_mapping][star][memory_bam_sort]} \
        -@ {snakemake.config[step_config][ngs_mapping][star][num_threads_bam_sort]} \
        -o $out
}}

# Merge and index STAR output
index_bam()
{{
    set -x

    out_bam=$1

    # Index resulting BAM file
    samtools index $out_bam

    # Build MD5 files
    pushd $(dirname $out_bam) && \
        md5sum $(basename $out_bam) > $(basename $out_bam).md5 && \
        md5sum $(basename $out_bam).bai > $(basename $out_bam).bai.md5 && \
    popd
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
        $(if [[ -n "{snakemake.config[step_config][ngs_mapping][star][quant_mode]}" ]]; then \
            echo --quantMode {snakemake.config[step_config][ngs_mapping][star][quant_mode]}
        fi) \
        --runThreadN {snakemake.config[step_config][ngs_mapping][star][num_threads_align]}

    >&2 ls -lhR $TMPDIR

    postproc_bam $TMPDIR/pre.d/out.Aligned.out.sam {snakemake.output.bam}

    if [[ "{snakemake.config[step_config][ngs_mapping][star][quant_mode]}" = *"TranscriptomeSAM"* ]]; then
        out_tx=$(echo "{snakemake.output.bam}" | sed -e "s/\\.bam$/\\.toTranscriptome\\.bam/")
        postproc_bam $TMPDIR/pre.d/out.Aligned.toTranscriptome.out.bam ${{out_tx}}
    fi
}}

# Perform Alignment -------------------------------------------------------------------------------

# Run STAR
run_star

index_bam {snakemake.output.bam}

if [[ "{snakemake.config[step_config][ngs_mapping][star][quant_mode]}" = *"TranscriptomeSAM"* ]]; then
    out_tx=$(echo "{snakemake.output.bam}" | sed -e "s/\\.bam$/\\.toTranscriptome\\.bam/")
    index_bam ${{out_tx}}
fi

if [[ "{snakemake.config[step_config][ngs_mapping][star][quant_mode]}" = *"GeneCounts"* ]]; then
    out_gc=$(echo "{snakemake.output.bam}" | sed -e "s/\\.bam$/\\.GeneCounts\\.tab/")
    mv $TMPDIR/pre.d/out.ReadsPerGene.out.tab ${{out_gc}}
    md5sum ${{out_gc}} > ${{out_gc}}.md5
fi

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

# Build MD5 files for the reports
md5sum {snakemake.output.report_bamstats_txt} > {snakemake.output.report_bamstats_txt_md5}
md5sum {snakemake.output.report_flagstats_txt} >{snakemake.output.report_flagstats_txt_md5}
md5sum {snakemake.output.report_idxstats_txt} > {snakemake.output.report_idxstats_txt_md5}

# Additional logging for transparency & reproducibility
if [[ ! -z "{snakemake.log.log}" ]]; then
    # Logging: Save a copy this wrapper (with the pickle details in the header)
    cp {this_file} $(dirname {snakemake.log.log})/wrapper_star.py
    # Logging: Save a permanent copy of the environment file used
    cp $(dirname {this_file})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper_star.yaml
    # Logging: parameters used by STAR
    cp $TMPDIR/pre.d/out.Log.out $(dirname {snakemake.log.log})/Log.out
    # Logging: STAR alignment summary
    cp $TMPDIR/pre.d/out.Log.final.out $(dirname {snakemake.log.log})/Log.final.out
fi

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
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

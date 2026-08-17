from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

out_gc = snakemake.output.get("gene_counts", "__dummy__")
out_sj = snakemake.output.get("junctions", "__dummy__")
out_tx = snakemake.output.get("transcriptome", "__dummy__")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = args["input"]["reads_left"]
reads_right = args["input"].get("reads_right", "")

ShellWrapper(snakemake).run(
    r"""
set -x

# Some More Preparation ---------------------------------------------------------------------------

mkdir -p $TMPDIR/tmp.d $TMPDIR/pre.d

# Define some global shortcuts
INDEX={args[path_index]}

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
        -m {args[memory_bam_sort]} \
        -@ {args[num_threads_bam_sort]} \
        -
}}

# Merge and index STAR output
index_bam()
{{
    set -x

    out_bam=$1

    # Index resulting BAM file
    samtools index $out_bam

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
    if [[ "{args[trim_adapters]}" == "True" ]]; then
        trim_cmd="\"trimadap-mt -p {args[num_threads_trimming]}\""
    else
        trim_cmd="zcat"
    fi

    quant_mode=""
    if [[ -n "{args[features]}" ]]
    then
        quant_mode="$quant_mode GeneCounts"
    fi
    if [[ "{args[transcriptome]}" = "True" ]]
    then
        quant_mode="$quant_mode TranscriptomeSAM"
    fi

    STAR \
        --readFilesIn ${{left_files}} ${{right_files}} \
        {args[raw_star_options]} \
        $rg_args \
        --readFilesCommand ${{trim_cmd}} \
        --alignIntronMax {args[align_intron_max]} \
        --alignIntronMin {args[align_intron_min]} \
        --alignMatesGapMax {args[align_mates_gap_max]} \
        --alignSJDBoverhangMin {args[align_sjdb_overhang_min]} \
        --alignSJoverhangMin {args[align_sj_overhang_min]} \
        --genomeDir {args[path_index]} \
        --genomeLoad {args[genome_load]} \
        --outFileNamePrefix $TMPDIR/pre.d/out. \
        --outFilterIntronMotifs {args[out_filter_intron_motifs]} \
        --outFilterMismatchNmax {args[out_filter_mismatch_n_max]} \
        --outFilterMismatchNoverLmax {args[out_filter_mismatch_n_over_l_max]} \
        --outFilterMultimapNmax {args[out_filter_multimap_n_max]} \
        --outFilterType {args[out_filter_type]} \
        --outSAMstrandField {args[out_sam_strand_field]} \
        --outSAMunmapped $(if [[ "{args[include_unmapped]}" == "True" ]]; then \
                echo "Within"; \
            else
                echo "None"; \
            fi) \
        $(if [[ -n "$quant_mode" ]]; then \
            echo "--quantMode $quant_mode"
        fi) \
        $(if [[ -n "{args[features]}" ]]; then \
            echo --sjdbGTFfile "{args[features]}"
        fi) \
        $(if [[ "{args[mask_duplicates]}" == "True" ]]; then \
            echo " --outStd SAM " ; \
        else
            echo " --outSAMtype BAM SortedByCoordinate "; \
        fi) \
        --runThreadN {args[num_threads_align]}

    >&2 ls -lhR $TMPDIR
}}

# Perform Alignment -------------------------------------------------------------------------------

# Run STAR
if [[ "{args[mask_duplicates]}" == "True" ]]; then
    run_star | mask_duplicates | sort_by_coord > {snakemake.output.bam}
else
    run_star
    mv $TMPDIR/pre.d/out.Aligned.sortedByCoord.out.bam {snakemake.output.bam}
fi

index_bam {snakemake.output.bam}

mv $TMPDIR/pre.d/out.ReadsPerGene.out.tab {out_gc}
mv $TMPDIR/pre.d/out.SJ.out.tab {out_sj}

# Optional output: mapping on transcriptome -------------------------------------------------------

if [[ "{args[transcriptome]}" = "True" ]]; then
    if [[ "{args[mask_duplicates]}" == "True" ]]; then
        samtools view -h -S $TMPDIR/pre.d/out.Aligned.toTranscriptome.out.bam | mask_duplicates | samtools view -h -b - > {out_tx}
    else
        mv $TMPDIR/pre.d/out.Aligned.toTranscriptome.out.bam {out_tx}
    fi
fi

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

"""
)

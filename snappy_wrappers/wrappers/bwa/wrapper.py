# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for BWA: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
input_left = args["input"]["reads_left"]
input_right = args["input"].get("reads_right", "")

path_bwa_index = args["path_index"]
trim_adapters = args["trim_adapters"]
num_threads_trimming = args["num_threads_trimming"]
mask_duplicates = args["mask_duplicates"]
num_threads_bam_view = args["num_threads_bam_view"]
memory_bam_sort = args["memory_bam_sort"]
num_threads_bam_sort = args["num_threads_bam_sort"]
num_threads_align = args["num_threads_align"]
split_as_secondary = args["split_as_secondary"]

ShellWrapper(snakemake).run(
    r"""
set -x
mkdir -p $TMPDIR/tmp.d

# Define some global shortcuts
INDEX={path_bwa_index}

# Define left and right reads as Bash arrays
declare -a reads_left=({input_left})
declare -a reads_right=({input_right})

# Check whether we have paired reads
if [[ ${{#reads_right[*]}} -eq 0 ]]; then
    paired=0
else
    paired=1
fi

# Check that we either have single-ended reds
if [[ $paired -eq 1 ]] && [[ ${{#reads_right[*]}} -ne ${{#reads_left[*]}} ]]; then
    >&2 echo "Number of right and left reads must be the same but was"
    >&2 echo "  left:  $reads_left"
    >&2 echo "  right: $reads_right"
    exit 1
fi

# Function Definitions ----------------------------------------------------------------------------

# Adapter trimming
trim_adapters()
{{
    set -x

    if [[ "{trim_adapters}" == "True" ]]; then
        trimadap-mt -p {num_threads_trimming}
    else
        cat  # TODO: can we somehow remove this?
    fi
}}

# Duplicate masking
mask_duplicates()
{{
    set -x

    if [[ "{mask_duplicates}" == "True" ]]; then
        samblaster --addMateTags
    else
        cat  # TODO: can we somehow remove this?
    fi
}}

# Alignment postprocessing
postproc_bam()
{{
    set -x
    out=$1

    mask_duplicates \
    | samtools view \
        -u \
        -Sb \
        -@ {num_threads_bam_view} \
    | samtools sort \
        -T $TMPDIR/sort_bam \
        -m {memory_bam_sort} \
        -@ {num_threads_bam_sort} \
        -O BAM \
        -o $out
}}

# Function for running BWA-ALN
run_bwa_aln()
{{
    set -euo pipefail
    set -x

    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
        # Compute suffix array indices for BWA-ALN
        bwa aln -t {num_threads_align} $INDEX ${{reads_left[$i]}} >$TMPDIR/left.sai
        if [[ $paired -eq 1 ]]; then
            sai_right=$TMPDIR/right.sai
            fastq_right=${{reads_right[$i]}}

            bwa aln -t {num_threads_align} $INDEX $fastq_right >$sai_right
            bwa_cmd=sampe
        else
            sai_right=
            fastq_right=
            bwa_cmd=samse
        fi

        add_rg()
        {{
            if [[ "{snakemake.params.args[sample_name]}" != "" ]]; then
                samtools addreplacerg \
                    --input-fmt SAM \
                    --output-fmt SAM \
                    -r "@RG\tID:{snakemake.params.args[sample_name]}.$i\tSM:{snakemake.params.args[sample_name]}\tPL:{snakemake.params.args[platform]}" \
                    /dev/stdin \
                    /dev/stdout
            else
                cat
            fi
        }}

        # Perform SAMSE/SAMPE
        bwa $bwa_cmd \
            $INDEX \
            $TMPDIR/left.sai \
            $sai_right \
            ${{reads_left[$i]}} \
            $fastq_right \
        | add_rg \
        | samtools view \
            -b \
            -@ {num_threads_bam_view} \
            -o $TMPDIR/tmp.d/out.$i.bam
    done
}}

# Function for running BWA-MEM
run_bwa_mem()
{{
    set -x

    # Decide whether to write split reads as supplementary or secondary (-M means secondary)
    split_as_supp_flag=
    if [[ "{split_as_secondary}" == "True" ]]; then
        split_as_supp_flag="-M"
    fi

    for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
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

        if [[ ! -z "$fastq_right" ]]; then
            seqtk mergepe ${{reads_left[$i]}} $fastq_right
        else
            zcat ${{reads_left[$i] }}
        fi \
        | trim_adapters \
        | bwa mem \
            $INDEX \
            $split_as_supp_flag \
            $rg_arg \
            -p \
            -t {num_threads_align} \
            /dev/stdin \
        | samtools view \
            -b \
            -@ {num_threads_bam_view} \
            -o $TMPDIR/tmp.d/out.$i.bam
    done
}}

# Perform Alignment -------------------------------------------------------------------------------

# estimate read length from first 100k reads
avg_len=$({{ zcat --force -- ${{reads_left[0]}} || true; }} \
          | head -n 400000 \
          | awk '(NR % 4 == 2) {{ count += 1; totLen += length($0) }}
                 END {{ print int(totLen/count) }}')

# Switch to BWA-SAMPE/SAMSE for shorter reads and BWA-MEM for longer ones.
if [[ $avg_len -le 75 ]]; then
    run_bwa_aln
else
    run_bwa_mem
fi

# Move over a single output file but merge multiple ones
if [[ ${{#reads_left[@]}} -eq 1 ]]; then
    samtools view \
        -@ {num_threads_bam_view} \
        -h $TMPDIR/tmp.d/out.0.bam \
    | postproc_bam {snakemake.output.bam}
else
    # Create merged header
    for f in $TMPDIR/tmp.d/out.*.bam; do
        samtools view \
            -@ {num_threads_bam_view} \
            -H $f >${{f%.bam}}.hdr.sam
    done
    samtools merge $TMPDIR/merged.hdr.bam $TMPDIR/tmp.d/out.*.hdr.sam

    # Concatenate files
    samtools cat \
        -h $TMPDIR/merged.hdr.bam \
        $TMPDIR/tmp.d/out.*.bam \
    | samtools view -h \
    | postproc_bam {snakemake.output.bam}
fi

# Index resulting BAM file
samtools index {snakemake.output.bam}

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}
"""
)

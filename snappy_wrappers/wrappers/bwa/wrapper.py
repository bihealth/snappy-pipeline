# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for BWA: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
input_left = snakemake.params.args["input"]["reads_left"]
input_right = snakemake.params.args["input"].get("reads_right", "")

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
mkdir -p $TMPDIR/tmp.d

# Define some global shortcuts
INDEX={snakemake.config[step_config][ngs_mapping][bwa][path_index]}

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

    if [[ "{snakemake.config[step_config][ngs_mapping][bwa][trim_adapters]}" == "True" ]]; then
        trimadap-mt -p {snakemake.config[step_config][ngs_mapping][bwa][num_threads_trimming]}
    else
        cat  # TODO: can we somehow remove this?
    fi
}}

# Duplicate masking
mask_duplicates()
{{
    set -x

    if [[ "{snakemake.config[step_config][ngs_mapping][bwa][mask_duplicates]}" == "True" ]]; then
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
        -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_view]} \
    | samtools sort \
        -T $TMPDIR/sort_bam \
        -m {snakemake.config[step_config][ngs_mapping][bwa][memory_bam_sort]} \
        -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_sort]} \
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
        bwa aln -t {snakemake.config[step_config][ngs_mapping][bwa][num_threads_align]} $INDEX ${{reads_left[$i]}} >$TMPDIR/left.sai
        if [[ $paired -eq 1 ]]; then
            sai_right=$TMPDIR/right.sai
            fastq_right=${{reads_right[$i]}}

            bwa aln -t {snakemake.config[step_config][ngs_mapping][bwa][num_threads_align]} $INDEX $fastq_right >$sai_right
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
            -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_view]} \
            -o $TMPDIR/tmp.d/out.$i.bam
    done
}}

# Function for running BWA-MEM
run_bwa_mem()
{{
    set -x

    # Decide whether to write split reads as supplementary or secondary (-M means secondary)
    split_as_supp_flag=
    if [[ "{snakemake.config[step_config][ngs_mapping][bwa][split_as_secondary]}" == "True" ]]; then
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
            -t {snakemake.config[step_config][ngs_mapping][bwa][num_threads_align]} \
            /dev/stdin \
        | samtools view \
            -b \
            -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_view]} \
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
        -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_view]} \
        -h $TMPDIR/tmp.d/out.0.bam \
    | postproc_bam {snakemake.output.bam}
else
    # Create merged header
    for f in $TMPDIR/tmp.d/out.*.bam; do
        samtools view \
            -@ {snakemake.config[step_config][ngs_mapping][bwa][num_threads_bam_view]} \
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

# Build MD5 files
pushd $(dirname {snakemake.output.bam})
md5sum $(basename {snakemake.output.bam}) > $(basename {snakemake.output.bam}).md5
md5sum $(basename {snakemake.output.bam_bai}) > $(basename {snakemake.output.bam_bai}).md5
popd

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

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=$( echo ${{dst}} | sed '0,/output\//{{s/output\//work\//}}' )
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

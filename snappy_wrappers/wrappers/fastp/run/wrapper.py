# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for fastp: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
library_name = snakemake.params.args["library_name"]
input_left = snakemake.params.args["input"]["reads_left"]
input_right = (
    snakemake.params.args["input"]["reads_right"]
    if snakemake.params.args["input"]["reads_right"]
    else {}
)

input_path_left = list(input_left.keys())
prefix_left = [x["relative_path"] for x in input_left.values()]
filename_left = [x["filename"] for x in input_left.values()]

assert len(input_path_left) == len(prefix_left)
assert len(input_path_left) == len(filename_left)

if input_right:
    input_path_right = list(input_right.keys())
    prefix_right = [x["relative_path"] for x in input_right.values()]
    filename_right = [x["filename"] for x in input_right.values()]

    assert len(input_path_left) == len(input_path_right)
    assert len(prefix_left) == len(prefix_right)
    assert len(filename_left) == len(filename_right)
else:
    input_path_right = []
    prefix_right = []
    filename_right = []

this_file = __file__

this_step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][this_step]["fastp"]

shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper_fastp.py

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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
mkdir -p $TMPDIR/out $TMPDIR/report

# Define left and right reads as Bash arrays
declare -a reads_left=({input_path_left})
declare -a reads_right=({input_path_right})
declare -a prefix_left=({prefix_left})
declare -a prefix_right=({prefix_right})
declare -a filename_left=({filename_left})
declare -a filename_right=({filename_right})

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

# Start fastp script

for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
    in=${{reads_left[$i]}}
    out=$TMPDIR/out/${{prefix_left[$i]}}/${{filename_left[$i]}}
    unpaired=$TMPDIR/rejected/${{prefix_left[$i]}}/${{filename_left[$i]}}
    failed=$TMPDIR/rejected/${{prefix_left[$i]}}/failed_${{filename_left[$i]}}
    mkdir -p $(dirname $out) $(dirname $unpaired)

    if [[ $paired -eq 1 ]]; then
        in2=${{reads_right[$i]}}
        out2=$TMPDIR/out/${{prefix_right[$i]}}/${{filename_right[$i]}}
        unpaired2=$TMPDIR/rejected/${{prefix_right[$i]}}/${{filename_right[$i]}}
        mkdir -p $(dirname $out2) $(dirname $unpaired2)
    fi

    # Find the base of filename (i.e. without fastq.gz (or derivatives) extension)
    report="${{prefix_left[$i]}}/${{filename_left[$i]}}"
    if [[ $paired -eq 1 ]]; then
        report="${{report}}_${{filename_right[$i]}}"
    fi
    mkdir -p $TMPDIR/report/$report

    fastp \
        --in1 $in --out1 $out \
        $(if [[ $paired -eq 1 ]]; then \
            echo --unpaired1 $unpaired --in2 $in2 --out2 $out2 --unpaired2 $unpaired2
        fi) \
        --failed_out $failed \
        --json $TMPDIR/report/$report/report.json \
        --html $TMPDIR/report/$report/report.html \
        --thread {config[num_threads]} \
        --trim_front1 {config[trim_front1]}                                 \
        --trim_tail1 {config[trim_tail1]}                                   \
        --max_len1 {config[max_len1]}                                       \
        --trim_front2 {config[trim_front2]}                                 \
        --trim_tail2 {config[trim_tail2]}                                   \
        --max_len2 {config[max_len2]}                                       \
        $(if [[ "{config[dedup]}" = "True" ]]; then \
            echo --dedup --dup_calc_accuracy {config[dup_calc_accuracy]}
        else \
            echo --dont_eval_duplication
        fi) \
        $(if [[ "{config[trim_poly_g]}" = "True" ]]; then \
            echo --trim_poly_g --poly_g_min_len {config[poly_g_min_len]}
        else \
            echo --disable_trim_poly_g
        fi) \
        $(if [[ "{config[trim_poly_x]}" = "True" ]]; then \
            echo --trim_poly_x --poly_x_min_len {config[poly_x_min_len]}
        fi) \
        $(if [[ "{config[cut_front]}" = "True" ]]; then \
            echo --cut_front --cut_front_window_size {config[cut_front_window_size]} --cut_front_mean_quality {config[cut_front_mean_quality]}
        fi) \
        $(if [[ "{config[cut_tail]}" = "True" ]]; then \
            echo --cut_tail --cut_tail_window_size {config[cut_tail_window_size]} --cut_tail_mean_quality {config[cut_tail_mean_quality]}
        fi) \
        $(if [[ "{config[cut_right]}" = "True" ]]; then \
            echo --cut_right --cut_right_window_size {config[cut_right_window_size]} --cut_right_mean_quality {config[cut_right_mean_quality]}
        fi) \
        $(if [[ "{config[disable_quality_filtering]}" = "False" ]]; then \
            echo --qualified_quality_phred {config[qualified_quality_phred]}
            echo --unqualified_percent_limit {config[unqualified_percent_limit]}
            echo --n_base_limit {config[n_base_limit]}
            echo --average_qual {config[average_qual]}
        else \
            echo --disable_quality_filtering
        fi) \
        $(if [[ "{config[disable_length_filtering]}" = "False" ]]; then \
            echo --length_required {config[length_required]}
            echo --length_limit {config[length_limit]}
        else \
            echo --disable_length_filtering
        fi) \
        $(if [[ "{config[low_complexity_filter]}" = "True" ]]; then \
            echo --low_complexity_filter --complexity_threshold {config[complexity_threshold]}
        fi) \
        $(if [[ -n "{config[filter_by_index1]}" ]]; then \
            echo --filter_by_index1 {config[filter_by_index1]} --filter_by_index_threshold {config[filter_by_index_threshold]}
        fi) \
        $(if [[ -n "{config[filter_by_index2]}" ]]; then \
            echo --filter_by_index2 {config[filter_by_index2]} --filter_by_index_threshold {config[filter_by_index_threshold]}
        fi) \
        $(if [[ "{config[correction]}" = "True" ]]; then \
            echo --correction
            echo --overlap_len_require {config[overlap_len_require]}
            echo --overlap_diff_limit {config[overlap_diff_limit]}
            echo --overlap_diff_percent_limit {config[overlap_diff_percent_limit]}
        fi) \
        $(if [[ "{config[umi]}" = "True" ]]; then \
            echo --umi
            echo --umi_loc {config[umi_loc]}
            echo --umi_len {config[umi_len]}
            echo --umi_prefix {config[umi_prefix]}
            echo --umi_skip {config[umi_skip]}
        fi) \
        $(if [[ "{config[overrepresentation_analysis]}" = "True" ]]; then \
            echo --overrepresentation_analysis
        fi)

    fns="$out $failed"
    if [[ $paired -eq 1 ]]; then
        fns="$unpaired $fns $out2 $unpaired2"
    fi
    fns=$(echo "$fns" | tr ' ' '\n')
    for fn in $fns ; do
        pushd $(dirname $fn)
        md5sum $fn > $fn.md5
        popd
    done

    pushd $TMPDIR/report/$report
    fns=$(ls)
    for fn in $fns ; do
        md5sum $fn > $fn.md5
    done
    popd

done

d=$(dirname {snakemake.output.out_done})
mkdir -p $d
cp -r $TMPDIR/out/* $d

d=$(dirname {snakemake.output.report_done})
mkdir -p $d
cp -r $TMPDIR/report/* $d

d=$(dirname {snakemake.output.rejected_done})
mkdir -p $d
cp -r $TMPDIR/rejected/* $d

touch {snakemake.output.out_done}
touch {snakemake.output.report_done}
touch {snakemake.output.rejected_done}

"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bbduk: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Eric Blanc <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
input_left = snakemake.params.args["input"]["reads_left"]
input_right = snakemake.params.args["input"].get("reads_right", "")

this_file = __file__

this_step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][this_step]["bbduk"]

shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper_bwa.py

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

# Start bbduk script

for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
    # Find the base of filename (i.e. without fastq.gz (or derivatives) extension)
    base=$(basename ${{reads_left[$i]}} | sed -e "s/\.f(ast)q(\.gz)?$//")
    if [[ $paired -eq 1 ]]; then
        # Common name between left & right reads (https://stackoverflow.com/questions/6973088/longest-common-prefix-of-two-strings-in-bash)
        base=$(printf "%s\n%s\n" $(basename ${{reads_left[$i]}}) $(basename ${{reads_right[$i]}}) | sed -e 'N;s/^\(.*\).*\n\1.*$/\1/' | sed -e "s/_R$//")
    fi

    in=${{reads_left[$i]}}
    out=$TMPDIR/out/$(basename $in)
    outm=$TMPDIR/out/rejected_$(basename $in)

    if [[ $paired -eq 1 ]]; then
        in2=${{reads_right[$i]}}
        out2=$TMPDIR/out/$(basename $in2)
        outm2=$TMPDIR/out/rejected_$(basename $in2)
    fi

    bbduk.sh \
        in=$in out=$out outm=$outm \
        $(if [[ $paired -eq 1 ]]; then \
            echo in2=$in2 out2=$out2 outm2=$outm2
        fi) \
        ref={config[adapter_sequences]} \
        stats=$TMPDIR/report/stats_$base.tsv \
        refstats=$TMPDIR/report/refstats_$base.tsv \
        rpkm=$TMPDIR/report/rpkm_$base.tsv \
        t={config[num_threads]} \
        {config[parameters]}

done

pushd $TMPDIR/out
fns=$(ls)
for fn in $fns; do
    md5sum $fn > $fn.md5
done
popd

pushd $TMPDIR/report
fns=$(ls)
for fn in $fns; do
    md5sum $fn > $fn.md5
done
popd

d=$(dirname {snakemake.output.out_done})
mkdir -p $d
cp $TMPDIR/out/* $d

d=$(dirname {snakemake.output.report_done})
mkdir -p $d
cp $TMPDIR/report/* $d

touch {snakemake.output.out_done}
touch {snakemake.output.report_done}

"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

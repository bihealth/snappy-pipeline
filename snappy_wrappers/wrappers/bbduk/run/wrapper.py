# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bbduk: Snakemake wrapper.py"""

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
config = snakemake.config["step_config"][this_step]["bbduk"]

shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper_bbduk.py

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
mkdir -p $TMPDIR/out $TMPDIR/rejected $TMPDIR/report

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

# Start bbduk script

for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
    in=${{reads_left[$i]}}
    out=$TMPDIR/out/${{prefix_left[$i]}}/${{filename_left[$i]}}
    outm=$TMPDIR/rejected/${{prefix_left[$i]}}/${{filename_left[$i]}}
    mkdir -p $(dirname $out) $(dirname $outm)

    if [[ $paired -eq 1 ]]; then
        in2=${{reads_right[$i]}}
        out2=$TMPDIR/out/${{prefix_right[$i]}}/${{filename_right[$i]}}
        outm2=$TMPDIR/rejected/${{prefix_right[$i]}}/${{filename_right[$i]}}
        mkdir -p $(dirname $out2) $(dirname $outm2)
    fi

    # Find the base of filename (i.e. without fastq.gz (or derivatives) extension)
    report="${{prefix_left[$i]}}/${{filename_left[$i]}}"
    if [[ $paired -eq 1 ]]; then
        report="${{report}}_${{filename_right[$i]}}"
    fi
    mkdir -p $TMPDIR/report/$report

    # Create dummy files for debugging
    # touch $out $outm $out2 $outm2
    # for fn in $(echo "stats refstats rpkm bhist qhist qchist aqhist bqhist lhist phist enthist" | tr ' ' '\n'); do
    #     touch $TMPDIR/report/$report/$fn.tsv
    # done

    bbduk.sh \
        in=$in out=$out outm=$outm \
        $(if [[ $paired -eq 1 ]]; then \
            echo in2=$in2 out2=$out2 outm2=$outm2
        fi) \
        ref=$(echo "{config[adapter_sequences]}" | tr ' ' ',') \
        stats=$TMPDIR/report/$report/stats.tsv       \
        refstats=$TMPDIR/report/$report/refstats.tsv \
        rpkm=$TMPDIR/report/$report/rpkm.tsv         \
        bhist=$TMPDIR/report/$report/bhist.tsv       \
        qhist=$TMPDIR/report/$report/qhist.tsv       \
        qchist=$TMPDIR/report/$report/qchist.tsv     \
        aqhist=$TMPDIR/report/$report/aqhist.tsv     \
        bqhist=$TMPDIR/report/$report/bqhist.tsv     \
        lhist=$TMPDIR/report/$report/lhist.tsv       \
        phist=$TMPDIR/report/$report/phist.tsv       \
        enthist=$TMPDIR/report/$report/enthist.tsv   \
        t={config[num_threads]} \
        interleaved={config[interleaved]}             \
        qin={config[qin]}                             \
        copyundefined={config[copyundefined]}         \
        nzo={config[nzo]}                             \
        qout={config[qout]}                           \
        statscolumns={config[statscolumns]}           \
        rename={config[rename]}                       \
        refnames={config[refnames]}                   \
        trd={config[trd]}                             \
        ordered={config[ordered]}                     \
        gcbins={config[gcbins]}                       \
        maxhistlen={config[maxhistlen]}               \
        histbefore={config[histbefore]}               \
        idbins={config[idbins]}                       \
        k={config[k]}                                 \
        rcomp={config[rcomp]}                         \
        maskmiddle={config[maskmiddle]}               \
        minkmerhits={config[minkmerhits]}             \
        minkmerfraction={config[minkmerfraction]}     \
        mincovfraction={config[mincovfraction]}       \
        hammingdistance={config[hammingdistance]}     \
        qhdist={config[qhdist]}                       \
        editdistance={config[editdistance]}           \
        hammingdistance2={config[hammingdistance2]}   \
        qhdist2={config[qhdist2]}                     \
        editdistance2={config[editdistance2]}         \
        forbidn={config[forbidn]}                     \
        removeifeitherbad={config[removeifeitherbad]} \
        trimfailures={config[trimfailures]}           \
        findbestmatch={config[findbestmatch]}         \
        skipr1={config[skipr1]}                       \
        skipr2={config[skipr2]}                       \
        ecco={config[ecco]}                           \
        ktrim={config[ktrim]}                         \
        $(if [[ -n "{config[kmask]}" ]] ; then \
            echo kmask={config[kmask]} 
        fi) \
        maskfullycovered={config[maskfullycovered]}   \
        ksplit={config[ksplit]}                       \
        mink={config[mink]}                           \
        qtrim={config[qtrim]}                         \
        trimq={config[trimq]}                         \
        minlength={config[minlength]}                 \
        mlf={config[mlf]}                             \
        minavgquality={config[minavgquality]}         \
        maqb={config[maqb]}                           \
        minbasequality={config[minbasequality]}       \
        maxns={config[maxns]}                         \
        mcb={config[mcb]}                             \
        ottm={config[ottm]}                           \
        tp={config[tp]}                               \
        tbo={config[tbo]}                             \
        strictoverlap={config[strictoverlap]}         \
        minoverlap={config[minoverlap]}               \
        mininsert={config[mininsert]}                 \
        tpe={config[tpe]}                             \
        forcetrimleft={config[forcetrimleft]}         \
        forcetrimright={config[forcetrimright]}       \
        forcetrimright2={config[forcetrimright2]}     \
        forcetrimmod={config[forcetrimmod]}           \
        restrictleft={config[restrictleft]}           \
        restrictright={config[restrictright]}         \
        mingc={config[mingc]}                         \
        maxgc={config[maxgc]}                         \
        tossjunk={config[tossjunk]}                   \
        swift={config[swift]}                         \
        chastityfilter={config[chastityfilter]}       \
        barcodefilter={config[barcodefilter]}         \
        $(if [[ -n "{config[barcodes]}" ]]; then \
            echo barcodes={config[barcodes]}  
        fi) \
        xmin={config[xmin]}                           \
        ymin={config[ymin]}                           \
        xmax={config[xmax]}                           \
        ymax={config[ymax]}                           \
        trimpolya={config[trimpolya]}                 \
        trimpolygleft={config[trimpolygleft]}         \
        trimpolygright={config[trimpolygright]}       \
        trimpolyg={config[trimpolyg]}                 \
        filterpolyg={config[filterpolyg]}             \
        entropy={config[entropy]}                     \
        entropywindow={config[entropywindow]}         \
        entropyk={config[entropyk]}                   \
        minbasefrequency={config[minbasefrequency]}   \
        entropytrim={config[entropytrim]}             \
        entropymask={config[entropymask]}             \
        entropymark={config[entropymark]}             \
        cardinality={config[cardinality]}             \
        cardinalityout={config[cardinalityout]}       \
        loglogk={config[loglogk]}                     \
        loglogbuckets={config[loglogbuckets]} 

    fns="$out $outm"
    if [[ $paired -eq 1 ]]; then
        fns="$fns $out2 $outm2"
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

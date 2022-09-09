# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GATK HaplotypeCaller: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

if snakemake.config["step_config"]["ngs_mapping"]["gatk_post_bam"]["do_realignment"]:
    bam_realigned = snakemake.output.bam_realigned
else:
    bam_realigned = ""
if snakemake.config["step_config"]["ngs_mapping"]["gatk_post_bam"]["do_recalibration"]:
    bam_recalibrated = snakemake.output.bam_recalibrated
else:
    bam_recalibrated = ""


shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

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

# Setup auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Get some shortcuts
realigned_infix={snakemake.config[step_config][ngs_mapping][gatk_post_bam][realigned_infix]}
recalibrated_infix={snakemake.config[step_config][ngs_mapping][gatk_post_bam][recalibrated_infix]}

# Build realigner target intervals file for all paths given as arguments
realigner_target_creator()
{{
    set -ex
    mkdir -p $TMPDIR/RealignerTargetCreator

    gatk_nonfree -Xmx40g -Djava.io.tmpdir=$TMPDIR \
        -T RealignerTargetCreator \
        -R {snakemake.config[static_data_config][reference][path]} \
        $(for path in {snakemake.config[step_config][ngs_mapping][gatk_post_bam][paths_known_sites]}; do \
            echo --known $path; \
        done) \
        $(for path in $*; do \
            echo -I $path; \
        done) \
        -o $TMPDIR/RealignerTargetCreator/forIndelRealigner.intervals
}}

# Perform indel realignment
indel_realigner()
{{
    set -ex
    mkdir -p $TMPDIR/IndelRealigner

    # Build the ".map" file
    for path in $*; do
        out_fname=$(basename $path .bam).$realigned_infix.bam
        echo -e "$(basename $path)\t$TMPDIR/IndelRealigner/$out_fname" \
        >> $TMPDIR/IndelRealigner/nWayOut.map
    done
    >&2 cat $TMPDIR/IndelRealigner/nWayOut.map

    # Perform indel realignment
    gatk_nonfree -Xmx50g -Djava.io.tmpdir=$TMPDIR \
        -T IndelRealigner \
        -R {snakemake.config[static_data_config][reference][path]} \
        $(for path in {snakemake.config[step_config][ngs_mapping][gatk_post_bam][paths_known_sites]}; do \
            echo --knownAlleles $path; \
        done) \
        $(for path in $*; do \
            echo -I $path; \
        done) \
        -targetIntervals $TMPDIR/RealignerTargetCreator/forIndelRealigner.intervals \
        --nWayOut $TMPDIR/IndelRealigner/nWayOut.map

    # Generate check sum files
    pushd $TMPDIR/IndelRealigner && \
        for fname in *.bai; do \
            mv $fname ${{fname%.bai}}.bam.bai; \
        done && \
        for fname in *.bam *.bai; do \
            md5sum $fname >$fname.md5; \
        done && \
        popd
}}

# Perform base recalibration
base_recalibrator()
{{
    set -ex
    mkdir -p $TMPDIR/BaseRecalibrator

    recalibrated_infix={snakemake.config[step_config][ngs_mapping][gatk_post_bam][recalibrated_infix]}
    for path in $*; do
        table_fname=$(basename $path .bam).$recalibrated_infix.grp
        out_fname=$(basename $path .bam).$recalibrated_infix.bam

        gatk_nonfree -Xmx50g -Djava.io.tmpdir=$TMPDIR \
            -T BaseRecalibrator \
            $(for path in {snakemake.config[step_config][ngs_mapping][gatk_post_bam][paths_known_sites]}; do \
                echo --knownSites $path; \
            done) \
            -R {snakemake.config[static_data_config][reference][path]} \
            -I $1 \
            -o $TMPDIR/BaseRecalibrator/$table_fname

        gatk_nonfree -Xmx50g -Djava.io.tmpdir=$TMPDIR \
            -T PrintReads \
            -R {snakemake.config[static_data_config][reference][path]} \
            -I $1 \
            -BQSR $TMPDIR/BaseRecalibrator/$table_fname \
            -o $TMPDIR/BaseRecalibrator/$out_fname
    done

    # Generate check sum files
    pushd $TMPDIR/BaseRecalibrator && \
        for fname in *.bai; do \
            mv $fname ${{fname%.bai}}.bam.bai; \
        done && \
        for fname in *.bam *.bai; do \
            md5sum $fname >$fname.md5; \
        done && \
        popd
}}

# Move realigned files to output directory, following the snakemake.input/snakemake.output paths
# and thus following the logic of the tool.
move_files()
{{
    set -ex
    result_dir=$1
    shift

    for path in $*; do
        output_dir=$(dirname $path)
        output_basename=$(basename $path .bam)
        mv $result_dir/$output_basename.* $output_dir
    done
}}

# Perform realignment if configured to do so
input="{snakemake.input}"
if [[ "{snakemake.config[step_config][ngs_mapping][gatk_post_bam][do_realignment]}" == "True" ]]; then
    realigner_target_creator $input
    indel_realigner $input
    output=$TMPDIR/IndelRealigner/*.bam
    input="$output"
fi

# Perform recalibration if configured to do so
if [[ "{snakemake.config[step_config][ngs_mapping][gatk_post_bam][do_recalibration]}" == "True" ]]; then
    base_recalibrator $input
fi

# Move step results files to output directory
if [[ "{snakemake.config[step_config][ngs_mapping][gatk_post_bam][do_realignment]}" == "True" ]]; then
    move_files $TMPDIR/IndelRealigner $(echo {bam_realigned} | tr ' ' ' \n' | grep '\.bam$')
fi
if [[ "{snakemake.config[step_config][ngs_mapping][gatk_post_bam][do_recalibration]}" == "True" ]]; then
    move_files $TMPDIR/BaseRecalibrator $(echo {bam_recalibrated} | tr ' ' ' \n' | grep '\.bam$')
fi
"""
)

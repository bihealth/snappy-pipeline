# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

import os

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

rscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "snappy-copywriter-run.R")
reference_folder = os.path.abspath(os.path.dirname(snakemake.input.gc))

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
# trap "rm -rf $TMPDIR" EXIT

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

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
cat <<"EOF" > $TMPDIR/run_copywriter_run.R
source("{rscript}")

runCopywriteR(
    normal.bam = "{snakemake.input.normal_bam}",
    tumor.bam = "{snakemake.input.tumor_bam}",
    destination.folder=".",
    binSize=as.numeric({args[bin_size]}),
    donorID="{snakemake.wildcards.library_name}",
    fullID="{snakemake.wildcards.library_name}",
    copywriter.params=list(
        reference.folder="{reference_folder}",
        capture.regions.file="{args[path_target_regions]}",
        workers=8
    )
)

warnings()
EOF

out_dir=$(readlink -f $(dirname {snakemake.output[0]}))
# and again
# need to remove "/out" because the R script writes to out and we would end up with out/out/
out_dir=$(dirname $out_dir)

cat $TMPDIR/run_copywriter_run.R

pushd $TMPDIR
mkdir -p out plots CNAprofiles
Rscript --vanilla run_copywriter_run.R

# -------------------------------------------------------------------------------------------------
# Move out output
#
cp -R * $out_dir
"""
)

shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
md5sum {snakemake.output.input} >{snakemake.output.input_md5}
md5sum {snakemake.output.segment} >{snakemake.output.segment_md5}
md5sum {snakemake.output.counts} >{snakemake.output.counts_md5}
md5sum {snakemake.output.log2} >{snakemake.output.log2_md5}
"""
)

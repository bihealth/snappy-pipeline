# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for converting exome baits & targets from bed to interval lists"""

import os
import re

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = args["reference"]
reference = re.sub("\.fa(sta)?(\.b?gz)?$", ".dict", reference)
assert os.path.exists(reference), "Missing dict of reference fasta"

baits = args["path_to_baits"]
targets = args.get("path_to_targets", "")

shell.executable("/bin/bash")

shell(
    r"""
set -x

d=$(ls $CONDA_PREFIX/share | grep picard)
picard_jar="$CONDA_PREFIX/share/$d/picard.jar"
if [[ ! -r $picard_jar ]]
then
    echo "Can't find picard jar"
    exit -1
fi

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

# Can't pipe to BedToIntervalList (https://github.com/broadinstitute/picard/issues/1890)
bed_to_interval_list() {{
    fn=$1
    f=$(basename $fn)
    if [[ $(od -x -N 2 $fn | head -n 1 | sed -e "s/.* //") = "8b1f" ]]
    then
        extract=zcat
    else
        extract=cat
    fi
    $extract $fn \
        | cut -f 1-3 \
        | bedtools sort -i - \
        | bedtools merge -i - \
        > $tmpdir/$f
    java -jar $picard_jar BedToIntervalList \
        -I $tmpdir/$f \
        -O /dev/stdout \
        -SD {reference}
}}

md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo $checksum
}}

bed_to_interval_list {baits} > {snakemake.output.baits}
md5 {snakemake.output.baits} > {snakemake.output.baits}.md5

if [[ -n "{targets}" ]]
then
    bed_to_interval_list {targets} > {snakemake.output.targets}
    md5 {snakemake.output.targets} > {snakemake.output.targets}.md5
else
    ln -rs {snakemake.output.baits} {snakemake.output.targets}
    ln -rs {snakemake.output.baits}.md5 {snakemake.output.targets}.md5
fi
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for converting exome baits & targets from bed to interval lists"""

import os
import re

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = args["reference"]
reference = re.sub(r"\.fa(sta)?(\.b?gz)?$", ".dict", reference)
assert os.path.exists(reference), "Missing dict of reference fasta"

baits = args["path_to_baits"]
targets = args.get("path_to_targets", "")

shell.executable("/bin/bash")

ShellWrapper(snakemake).run(
    r"""
set -x

d=$(ls $CONDA_PREFIX/share | grep picard)
picard_jar="$CONDA_PREFIX/share/$d/picard.jar"
if [[ ! -r $picard_jar ]]
then
    echo "Can't find picard jar"
    exit -1
fi

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


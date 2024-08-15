# -*- coding: utf-8 -*-
"""CUBI BEDVenn: Snakemake wrapper.py"""

import itertools
import textwrap

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def powerset(iterable):
    """powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"""
    s = list(iterable)
    return itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s) + 1))


all_ids = set(snakemake.input.keys())
id_sets = tuple(map(set, powerset(all_ids)))

PREFIX = r"""
set -x

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

export TMPDIR=$(mktemp -d)
"""

unique_cmd = "| sort -u" if snakemake.params["args"]["unique"] else ""

chunks = []

for id_set in id_sets:
    if not id_set:
        continue  # skip empty set
    labels = "\\t".join(sorted(id_set))
    intersect_files = [snakemake.input[name] for name in id_set]
    subtract_files = [
        snakemake.input[name] for name in snakemake.input.keys() if name not in id_set
    ]
    chunks.append(
        textwrap.dedent(
            r"""
    # Hack: get back bin directory of base/root environment.
    export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

    cut -f 1-3 {first} >$TMPDIR/input_first.txt
    if [[ -z "{intersect_files}" ]]; then
        cp $TMPDIR/input_first.txt $TMPDIR/intersected.txt
    else
        i=0
        cp $TMPDIR/input_first.txt $TMPDIR/intersected.0.txt
        for path in {intersect_files}; do
            let "i=$i+1"
            cut -f 1-3 $path > $TMPDIR/input_second.txt
            snappy-bed_jaccard_operations \
                --input-first $TMPDIR/intersected.$(($i - 1)).txt \
                --input-second $TMPDIR/input_second.txt \
                --threshold {{snakemake.params[args][jaccard_threshold]}} \
                --output-file /dev/stdout \
                --operation intersect \
            {unique_cmd} \
            > $TMPDIR/intersected.$i.txt
        done
        cp $TMPDIR/intersected.$i.txt $TMPDIR/intersected.txt
    fi

    if [[ -z "{subtract_files}" ]]; then
        cp $TMPDIR/intersected.txt $TMPDIR/subtracted.txt
    else
        cut -f 1-3 {subtract_files} >$TMPDIR/input_second.txt
        snappy-bed_jaccard_operations \
            --input-first $TMPDIR/intersected.txt \
            --input-second $TMPDIR/input_second.txt \
            --threshold {{snakemake.params[args][jaccard_threshold]}} \
            --output-file $TMPDIR/subtracted.txt \
            --operation subtract
    fi

    echo -e "{labels}\t$(wc -l $TMPDIR/subtracted.txt | awk '{{{{ print $1 }}}}')" \
    >>$TMPDIR/overlap.txt
    """
        ).format(
            first=intersect_files[0],
            intersect_files=" ".join(intersect_files[1:]),
            subtract_files=" ".join(subtract_files),
            unique_cmd=unique_cmd,
            labels=labels,
        )
    )

SUFFIX = r"""
export MPLBACKEND="agg"
snappy-quickvenn \
    --input-shared-counts $TMPDIR/overlap.txt \
    --output-image {snakemake.output.image}
pushd $(dirname {snakemake.output.image}) && \
    md5sum $(basename {snakemake.output.image}) >$(basename {snakemake.output.image}).md5 && \
    popd

cp $TMPDIR/overlap.txt {snakemake.output.counts}
pushd $(dirname {snakemake.output.counts}) && \
    md5sum $(basename {snakemake.output.counts}) >$(basename {snakemake.output.counts}).md5 && \
    popd
"""

shell.executable("/bin/bash")

shell(PREFIX + "\n".join(chunks) + SUFFIX)

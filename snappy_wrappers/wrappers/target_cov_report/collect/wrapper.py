# -*- coding: utf-8 -*-
"""CUBI Snakemake wrapper for collecting target coverage reports.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

OUT_REPORT={snakemake.output.txt}

first=1
for path in $(echo {snakemake.input} | tr ' ' '\n' | sort); do
    if [[ $first -eq 1 ]]; then
        first=0
        echo -en "sample\t" \
        > $OUT_REPORT

        egrep '^SN|^AMC' $path \
        | cut -f 2 \
        | tr '\n' '\t' \
        | perl -p -e 's/\t$/\n/' \
        | perl -p -e 's/(\d+)/% >=\1x/g' \
        >> $OUT_REPORT
    fi

    echo -en "$(basename $path .txt | sed -e 's/bwa.//')\t" \
    >> $OUT_REPORT

    egrep '^SN|^AMC' $path \
    | cut -f 3 \
    | tr '\n' '\t' \
    | sed -e 's/\t$/\n/' \
    >> $OUT_REPORT
done

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)

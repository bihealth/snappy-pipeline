# -*- coding: utf-8 -*-
"""CUBI Snakemake wrapper for collecting target coverage reports.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

OUT_REPORT={snakemake.output.txt}

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

md5sum {snakemake.output.txt} > {snakemake.output.txt_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
"""
)

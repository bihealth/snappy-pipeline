# -*- coding: utf-8 -*-
"""Snakemake wrapper for collecting the MSDN.
"""

import os.path
import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

assert len(snakemake.input.gatk_hc) == len(snakemake.input.gatk_ug)

with tempfile.TemporaryDirectory() as tmpdir:
    # Write paths to input files into temporary directory.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    path_gatk_hc = os.path.join(tmpdir, "gatk_hc")
    with open(path_gatk_hc, "wt") as tmpf:
        print("\n".join(snakemake.input.gatk_hc), file=tmpf)
    path_gatk_ug = os.path.join(tmpdir, "gatk_ug")
    with open(path_gatk_ug, "wt") as tmpf:
        print("\n".join(snakemake.input.gatk_ug), file=tmpf)

    shell(
        r"""
    set -x
    set -euo pipefail

    BP_COUNT=20

    TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    pairs=$(paste {path_gatk_hc} {path_gatk_ug})

    get_sample()
    {{
        txt=$1

        echo $(basename $txt .summary.txt) \
        | sed -e 's/.*de_novos_hard\.//g' \
        | cut -d - -f 1
    }}

    extract()
    {{
        txt=$1
        sample=$(get_sample $txt)

        tail -n +2 $txt \
        | cut -f 1,2,4,5 \
        | sed -e "s/^/$sample\t/"
    }}

    cluster()
    {{
        awk -v BP_COUNT=$BP_COUNT -F $'\t' '
            BEGIN {{
                OFS=FS;
                count = 0;
                prev = 0;
                prevIndex = 0;
                prevChrom = 0;
                prevPos = 0;
                printed = 0;
            }}

            {{
                $2 = "chr" $2;
                if (prevChrom == $2 && $3 - prevPos < BP_COUNT) {{
                    if (!printed) {{
                        count++;
                        print "--";
                        print $1 "-" count, prev;
                    }}
                    print $1 "-" count, $0;
                    printed = 1;
                }} else {{
                    printed = 0;
                }}

                prev = $0;
                prevIndex = $1;
                prevChrom = $2;
                prevPos = $3;
            }}
        '
    }}

    (
        echo -e "cluster\tindex\tchrom\tpos\tref\talt"

        for pair in $pairs; do
            f_hc=$(echo $pair | cut -d ';' -f 1)
            f_ug=$(echo $pair | cut -d ';' -f 2)
            sample=$(get_sample $f_hc)

            comm -1 -2 <(extract $f_hc | sort) <(extract $f_ug | sort) \
            | sort -k2,2 -k3,3n \
            | cluster \
            > $TMPDIR/lines.txt

            if [[ -s $TMPDIR/lines.txt ]]; then
                cat $TMPDIR/lines.txt
            fi
        done
    ) > {snakemake.output.txt}

    pushd $(dirname {snakemake.output.txt})
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
    """
    )

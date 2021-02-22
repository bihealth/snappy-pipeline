# -*- coding: utf-8 -*-
"""Wrapper for converting bcftools roh output to BED files.

One BED file is created for each sample in the .txt file but created in a "tmp" directory.
A later link-out step in the Snakefile will create the real output file.  This is because we
cannot have output functions.
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell(
    r"""
set -x

base=$(basename $(dirname $(dirname {snakemake.output.done})) | cut -d . -f 1-3)
mkdir -p $(dirname {snakemake.output.done})/tmp

# Convert the TXT output of bcftools roh to BED files ---------------------------------------------

zcat {snakemake.input.txt} \
| awk \
    -v base=$(dirname {snakemake.output})/tmp/$base. \
    -F $'\t' '
    BEGIN {{ OFS = FS }}
    /^RG/ {{
        print $3, $4 - 1, $5, "qual:" $8 "__markers:" $7 "__length:" $6 "bp" > base $2 ".bed";
    }}
    '

# bgzip and tabix-index the BED files -------------------------------------------------------------

pushd $(dirname {snakemake.output})/tmp

for f in *.bed; do
    bgzip $f
    tabix -f $f.gz

    md5sum $f.gz >$f.gz.md5
    md5sum $f.gz.tbi >$f.gz.tbi.md5
done
"""
)

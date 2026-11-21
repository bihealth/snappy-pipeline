# -*- coding: utf-8 -*-
"""Wrapper for running "popdel call" on the windows of one part of the genome.

Each line of the ``region`` BED file is a padded window; popdel calls each window on its own.
"""

import tempfile
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"


with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.input.profile), file=tmpf)
    tmpf.flush()
    # Actually run the script.
    ShellWrapper(snakemake).run(
        r"""
for name in $(cat {tmpf.name}); do
    echo $name >> $TMPDIR/profiles.txt
done

# BED windows (0-based, half-open) become popdel regions (1-based, closed)
i=0
while read -r chrom begin end; do
    i=$((i + 1))
    popdel call \
        -r "$chrom:$((begin + 1))-$end" \
        -o $TMPDIR/window.$i.vcf \
        {extra_arguments} \
        $TMPDIR/profiles.txt
    echo $TMPDIR/window.$i.vcf >> $TMPDIR/windows.txt
done < {snakemake.input.region}

bcftools concat --file-list $TMPDIR/windows.txt -o $TMPDIR/tmp.vcf

cat > $TMPDIR/header.txt << EOF
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
EOF

# Profile paths no longer carry a dotted prefix, so drop any directory part before splitting.
for sample in $(bcftools view --header-only $TMPDIR/tmp.vcf | grep '^#CHROM' | cut -f 10-); do
    echo -e "$sample\t$(basename $sample | rev | cut -d . -f 1 | rev)" >> $TMPDIR/samples.txt
done

bcftools annotate \
    --header-lines $TMPDIR/header.txt \
    $TMPDIR/tmp.vcf \
    | awk -F $'\t' \
        'BEGIN {{ OFS = FS }}
     /^#/ {{ print $0; }}
     /^[^#]/ {{ $8 = $8 ";SVMETHOD=POPDELv1.1.0"; print; }}' \
    | bcftools reheader \
        --samples $TMPDIR/samples.txt \
    | bcftools sort \
        -T $TMPDIR \
        -O z \
        -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}
"""
    )

# -*- coding: utf-8 -*-
"""Wrapper for running "popdel call"."""

import tempfile
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})


def unescape_dots_dashes(s: str) -> str:
    """Unescape dots and dashes from double-underscore constructs."""
    return s.replace("__hyphen__", "-").replace("__dot__", ".").replace("__under__", "_")


chrom = unescape_dots_dashes(args["chrom"])


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
    echo $name >>$TMPDIR/profiles.txt
done

popdel call \
    -r {chrom}:{args[begin]}-{args[end]} \
    -o $TMPDIR/tmp.vcf \
    $TMPDIR/profiles.txt

cat >$TMPDIR/header.txt <<EOF
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
EOF

for sample in $(bcftools view --header-only $TMPDIR/tmp.vcf | grep '^#CHROM' | cut -f 10-); do
    echo -e "$sample\t$(echo $sample | rev | cut -d . -f 1 | rev)" >>$TMPDIR/samples.txt
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

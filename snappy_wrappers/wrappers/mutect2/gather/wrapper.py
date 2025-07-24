# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2 gather."""

from snakemake import shell

__author__ = "Till Hartmann <till.hartmann@bih-charite.de>"

output = snakemake.output
input = snakemake.input

shell(r"""
set -x

# Concatenate raw calls vcfs & index result ----------------------
bcftools concat \
         --allow-overlaps \
         -d none \
            -o {output.raw} \
               -O z \
    {input.raw}
tabix -f {output.raw}

# Concatenate stats with GATK tool -------------------------------
stats=$(echo "{input.stats}" | sed -e "s/ / -stats /g")
gatk MergeMutectStats -stats $stats -O {output.stats}

# Create orientation model from all F1R2 files -------------------
orientation=$(echo "{input.f1r2}" | sed -e "s/ / -I /g")
gatk LearnReadOrientationModel -I $orientation -O {output.orientation}

# Compute md5 sums -----------------------------------------------
pushd $(dirname {output.raw})
for f in *; do
if [[ -f "$f" ]] && [[ $f != *.md5 ]]; then
md5sum $f > $f.md5;
fi;
done
popd
""")

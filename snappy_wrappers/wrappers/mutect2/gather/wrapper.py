# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2 gather."""

from snakemake import shell

__author__ = "Till Hartmann <till.hartmann@bih-charite.de>"

output = snakemake.output
input = snakemake.input

if getattr(input, "stats", []):
    assert getattr(output, "stats"), "Missing output file for stats"
    stats = "-stats ".join([str(stat) for stat in input.stats])
    output_stats = output.stats
else:
    stats = ""
    output_stats = ""

if getattr(input, "f1r2", []):
    assert getattr(output, "orientation"), "Missing output file for orientation"
    orientation = "-I ".join([str(f1r2) for f1r2 in input.f1r2])
    output_orientation = output.orientation
else:
    orientation = ""
    output_orientation = ""

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
if [[ -n "{stats}" ]]
then
    gatk MergeMutectStats -stats {stats} -O {output_stats}
fi

# Create orientation model from all F1R2 files -------------------
if [[ -n "{orientation}" ]]
then
    gatk LearnReadOrientationModel -I {orientation} -O {output_orientation}
fi

# Compute md5 sums -----------------------------------------------
pushd $(dirname {output.raw})
for f in *; do
if [[ -f "$f" ]] && [[ $f != *.md5 ]]; then
md5sum $f > $f.md5;
fi;
done
popd
""")

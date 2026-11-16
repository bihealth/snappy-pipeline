# -*- coding: utf-8 -*-
"""The VarFish table of the structural variants of a pedigree.

The callers' VCFs are prepared as before (Manta inversions, the SVLEN header of some callers,
the MELT header, a GT header, normalisation), then ``mehari annotate strucvars`` writes the table.
"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args

# Get shortcut to "fix_manta_invs.py" postprocessing script
fix_manta_invs = os.path.join(os.path.dirname(__file__), "fix_manta_invs.py")

ShellWrapper(snakemake).run(
    r"""
# Fix the Manta inversions
i=0
for vcf in {snakemake.input.vcf}; do
    let "i=$i+1"
    num=$(printf %03d $i)

    python3 {fix_manta_invs} \
        --reference-fasta {snakemake.input.reference} \
        --input-vcf $vcf \
        --output-vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf
    bcftools sort -o $TMPDIR/fixed_bnd_to_inv.$num.vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf

    # Fixup SVLEN=1 to SVLEN=.
    sed -i -e 's/ID=SVLEN,Number=1/ID=SVLEN,Number=./g' $TMPDIR/fixed_bnd_to_inv.$num.vcf
    # Fixup MELT header
    sed -i -e "s/seperated by '..'/separated by '\\\\\\\\|'/" $TMPDIR/fixed_bnd_to_inv.$num.vcf

    # Add the missing "GT" tag
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
        > $TMPDIR/header.gt.txt

    # Annotate and normalise the VCF files.
    bcftools annotate \
        -h $TMPDIR/header.gt.txt \
        $TMPDIR/fixed_bnd_to_inv.$num.vcf \
        | bcftools norm \
            -m -any \
            -c w \
            --fasta-ref {snakemake.input.reference} \
            -O z \
            -o $TMPDIR/final_for_import.$num.vcf.gz
    tabix -s1 -b2 -e2 -f $TMPDIR/final_for_import.$num.vcf.gz
done

mehari annotate strucvars \
    --assembly {args[assembly]} \
    --rng-seed {args[rng_seed]} \
    --path-input-ped {snakemake.input.ped} \
    $(for p in $TMPDIR/final_for_import.*.vcf.gz; do
        echo --path-input-vcf $p
    done) \
    {extra_arguments} \
    --output {snakemake.output.gts} \
    --output-format tsv
"""
)

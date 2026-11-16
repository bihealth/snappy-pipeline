# -*- coding: utf-8 -*-
"""FDR control of varlociraptor calls, once per variant type, concatenated into one VCF."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
events = " ".join(args["events"])
variant_types = " ".join(args["variant_types"])

ShellWrapper(snakemake).run(
    r"""
for variant_type in {variant_types}; do
    varlociraptor filter-calls control-fdr \
        {snakemake.input.vcf} \
        --mode {args[mode]} \
        --events {events} \
        --fdr {args[fdr]} \
        --var $variant_type \
        {extra_arguments} \
        > $TMPDIR/$variant_type.bcf
    bcftools index --force $TMPDIR/$variant_type.bcf
done

bcftools concat \
    --allow-overlaps \
    $TMPDIR/*.bcf \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type z \
        --output {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)

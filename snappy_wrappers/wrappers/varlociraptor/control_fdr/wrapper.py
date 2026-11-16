# -*- coding: utf-8 -*-
"""FDR control of varlociraptor calls, once per variant type, concatenated into one VCF.

In ``tag`` mode, the records that fail the control stay in the VCF with the task name in
``FILTER``: the passing records are a site annotation that marks the others with a flag, which
``bcftools filter`` turns into the soft filter before the flag is dropped again.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
events = " ".join(args["events"])
variant_types = " ".join(args["variant_types"])
mode = args["mode"]
filter_name = args.get("filter_name")  # the task name, only in tag mode

ShellWrapper(snakemake).run(
    r"""
for variant_type in {variant_types}; do
    varlociraptor filter-calls control-fdr \
        {snakemake.input.vcf} \
        --mode {args[fdr_mode]} \
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
        --output-type b \
        --output $TMPDIR/passing.bcf

if [[ "{mode}" == "tag" ]]; then
    # The passing records mark the others with a flag, which becomes the soft filter
    bcftools index --force $TMPDIR/passing.bcf
    bcftools annotate \
        --annotations $TMPDIR/passing.bcf \
        --mark-sites -FAILED_FDR \
        {snakemake.input.vcf} \
        | bcftools filter \
            --soft-filter {filter_name} \
            --mode + \
            --exclude 'INFO/FAILED_FDR=1' \
        | bcftools annotate \
            --remove INFO/FAILED_FDR \
            --output-type z \
            --output {snakemake.output.vcf}
else
    bcftools view \
        --output-type z \
        --output {snakemake.output.vcf} \
        $TMPDIR/passing.bcf
fi
tabix -f {snakemake.output.vcf}
"""
)

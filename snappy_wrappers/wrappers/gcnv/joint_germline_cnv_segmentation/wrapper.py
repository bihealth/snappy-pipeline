from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# Run actual tools --------------------------------------------------------------------------------

gatk JointGermlineCNVSegmentation \
    --reference {args[reference]} \
    $(for vcf in {snakemake.input.vcf}; do echo --variant $vcf; done) \
    --model-call-intervals {snakemake.input.interval_list} \
    --pedigree {snakemake.input.ped} \
    --output {snakemake.output.vcf}
"""
)

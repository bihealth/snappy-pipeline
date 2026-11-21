from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

paths_tsv = " ".join(snakemake.input.tsv)

ShellWrapper(snakemake).run(
    r"""
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

mkdir -p {snakemake.output.ploidy}

gatk DetermineGermlineContigPloidy \
    --model {snakemake.input.model}/ploidy-model \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --output {snakemake.output.ploidy} \
    --output-prefix ploidy

# The ploidy calls follow the PED files: X and Y their sex, the autosomes 2
python $(dirname {__file__})/../../../tools/gcnv_ploidy.py \
    --autosomes-diploid \
    {snakemake.output.ploidy} \
    {snakemake.input.ped}
"""
)

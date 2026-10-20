from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

ShellWrapper(snakemake).run(
    r"""
export TMPDIR=$(realpath $TMPDIR)
trap "rm -rf $TMPDIR" EXIT

jar=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)

java -jar $jar -T ReadBackedPhasing \
    -R {snakemake.input.reference} \
    -I {snakemake.input.bam} \
    -L {snakemake.input.vcf} \
    --variant {snakemake.input.vcf} \
    -o {snakemake.output.vcf}
"""
)

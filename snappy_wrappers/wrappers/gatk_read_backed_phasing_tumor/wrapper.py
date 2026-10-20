"""GATK 3 ReadBackedPhasing of a tumor library's combined somatic and germline variants with its
alignments, as pVACtools describes it for the proximal variants of neoepitope prediction"""

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
    -o $TMPDIR/phased.vcf

bgzip -c $TMPDIR/phased.vcf > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)

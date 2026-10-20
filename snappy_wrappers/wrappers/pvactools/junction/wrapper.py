import json

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

# The strandedness decision of the RNA library: 1 forward, 2 reverse
with open(snakemake.input.strandedness, "rt") as f:
    decision = str(json.load(f).get("decision"))
strandedness = {"1": "FR", "2": "RF"}.get(decision)
if strandedness is None:
    raise ValueError(f"regtools needs a stranded library, but the decision is {decision!r}")

ShellWrapper(snakemake).run(
    r"""
export TMPDIR=$(realpath $TMPDIR)
trap "rm -rf $TMPDIR" EXIT
vcf=$TMPDIR/normalized.vcf.gz

case "{snakemake.input.features}" in
    *.gz) zcat {snakemake.input.features} > $TMPDIR/features.gtf ;;
    *) ln -sr {snakemake.input.features} $TMPDIR/features.gtf ;;
esac

bcftools norm \
    --multiallelics -both \
    --output $vcf --output-type z --write-index=tbi \
    {snakemake.input.annotated}

regtools cis-splice-effects identify -s {strandedness} \
    -o {snakemake.output.junctions} \
    $vcf \
    {snakemake.input.bam} \
    {snakemake.input.reference} \
    $TMPDIR/features.gtf
"""
)

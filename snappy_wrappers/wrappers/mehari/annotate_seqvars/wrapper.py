# -*- coding: utf-8 -*-
"""The VarFish table of the small variants of a pedigree.

The calls are normalised (and cut to the exons, if given), annotated with ``mehari annotate
seqvars`` into a VCF in the job's temporary directory, and ``mehari postprocess varfish-seqvars``
writes the table from it with the pedigree and the HGNC table. The callers' INFO fields are
dropped first: the table takes none of them, and the strict VCF reader of the postprocessing
rejects header definitions of reserved keys that differ from the VCF specification, such as
varlociraptor's ``SVLEN`` and ``EVENT``. The pedigree file is a result, too, as ``varfish-cli``
imports it with the table.
"""

from typing import TYPE_CHECKING

from snappy_wrappers import utils
from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
config = args["config"]
assembly = config["assembly"]
mehari_options = " \\\n        ".join(utils.mehari_options(config))

transcript_dbs = snakemake.input.transcript_db
if isinstance(transcript_dbs, str):
    transcript_dbs = [transcript_dbs]
transcripts = " ".join("--transcripts {}".format(path) for path in transcript_dbs)

exon_bed = snakemake.input.get("exon_bed", "")
# Without exons, the whole VCF; bcftools takes "" as no region file only through a variable
regions = "--regions-file {}".format(exon_bed) if exon_bed else ""

ShellWrapper(snakemake).run(
    r"""
bcftools view \
    {regions} \
    {snakemake.input.vcf} \
    | bcftools norm \
        --multiallelics -any \
        --force \
        --fasta-ref {snakemake.input.reference} \
    | bcftools annotate \
        --remove INFO \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type z \
        --output $TMPDIR/normalized.vcf.gz
tabix -f $TMPDIR/normalized.vcf.gz

mehari annotate seqvars \
    {mehari_options} \
    {transcripts} \
    --frequencies {snakemake.input.frequency_db} \
    --clinvar {snakemake.input.clinvar_db} \
    {extra_arguments} \
    --reference {snakemake.input.reference} \
    --input $TMPDIR/normalized.vcf.gz \
    --output $TMPDIR/annotated.vcf.gz

mehari postprocess varfish-seqvars \
    --assembly {assembly} \
    --pedigree {snakemake.input.ped} \
    --hgnc {snakemake.input.hgnc_tsv} \
    --input $TMPDIR/annotated.vcf.gz \
    --output {snakemake.output.gts}

cp -H {snakemake.input.ped} {snakemake.output.ped}
"""
)

import os

from pathlib import Path
from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
path_exon_bed = args["path_exon_bed"]
reference = args["reference"]
transcript_db = args.get("transcript_db")
clinvar_db = args.get("clinvar_db")
frequency_db = args.get("frequency_db")
hgnc_tsv = args["hgnc_tsv"]

if not Path(transcript_db).exists(follow_symlinks=True):
    transcript_db = None
if not Path(clinvar_db).exists(follow_symlinks=True):
    clinvar_db = None
if not Path(frequency_db).exists(follow_symlinks=True):
    frequency_db = None
if not Path(hgnc_tsv).exists(follow_symlinks=True):
    raise ValueError(f"hgnc.tsv required for mehari tsv output but not found at {hgnc_tsv}")

if not (transcript_db or clinvar_db or frequency_db):
    raise ValueError(
        "At least one of the following databases must be provided: "
        "transcript_db, clinvar_db, frequency_db."
    )

transcript_db_param = f"--transcripts {transcript_db}" if transcript_db else ""
clinvar_db_param = f"--clinvar {clinvar_db}" if clinvar_db else ""
frequency_db_param = f"--frequencies {frequency_db}" if frequency_db else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# Extract around BED file, if given.  Otherwise, "just" normalize.
if [[ -n "{path_exon_bed}" ]] && [[ "{path_exon_bed}" != "None" ]]; then
    set -e
    bcftools view \
        -R {path_exon_bed} \
        {snakemake.input.vcf} \
    | bcftools norm \
        -m -any \
        --force \
        --fasta-ref {reference} \
    | bcftools sort -T $TMPDIR \
    | bgzip -c \
    > $TMPDIR/tmp.vcf.gz
    tabix -f $TMPDIR/tmp.vcf.gz
else
    set -e
    bcftools norm \
        -m -any \
        --force \
        --fasta-ref {reference} \
        {snakemake.input.vcf} \
    | bcftools sort -T $TMPDIR \
    | bgzip -c \
    > $TMPDIR/tmp.vcf.gz
    tabix -f $TMPDIR/tmp.vcf.gz
fi

# Perform Mehari sequence variant annotation.
# Note: The TSV export does NOT select MANE transcripts on its own,
# therefore --pick-transcript and --pick-transcript-mode are needed.
# This will NOT work properly when a tx-db with more than one source is used
mehari \
    annotate \
    seqvars \
    --reference {reference} \
    {transcript_db_param} {clinvar_db_param} {frequency_db_param} \
    --pick-transcript mane-select \
    --pick-transcript mane-plus-clinical \
    --pick-transcript length \
    --pick-transcript-mode first \
    --keep-intergenic \
    --path-input-ped {snakemake.input.ped} \
    --path-input-vcf $TMPDIR/tmp.vcf.gz \
    --hgnc {hgnc_tsv} \
    --path-output-tsv {snakemake.output.gts}

# FIXME this should really not be hardcoded.
cat <<EOF | gzip -c > {snakemake.output.db_infos}
genomebuild	db_name	release
GRCh37	clinvar	20210728
GRCh37	gnomad_exomes	r2.1.1
GRCh37	gnomad_genomes	r2.1.1
EOF

# Copy out PED file to output
cp -H {snakemake.input.ped} {snakemake.output.ped}
"""
)

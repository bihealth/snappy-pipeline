import os

from pathlib import Path
from snakemake.shell import shell

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

DEF_HELPER_FUNCS = r"""
compute-md5()
{
    if [[ $# -ne 2 ]]; then
        >&2 echo "Invalid number of arguments: $#"
        exit 1
    fi
    md5sum $1 \
    | awk '{ gsub(/.*\//, "", $2); print; }' \
    > $2
}
"""

shell(
    r"""
set -x

# Write files for reproducibility -----------------------------------------------------------------

{DEF_HELPER_FUNCS}

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
compute-md5 {snakemake.log.conda_list} {snakemake.log.conda_list_md5}
compute-md5 {snakemake.log.conda_info} {snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
compute-md5 {snakemake.log.wrapper} {snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
compute-md5 {snakemake.log.env_yaml} {snakemake.log.env_yaml_md5}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Run actual tools --------------------------------------------------------------------------------

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

# Compute MD5 sums on output files
compute-md5 {snakemake.output.db_infos} {snakemake.output.db_infos_md5}
compute-md5 {snakemake.output.gts} {snakemake.output.gts_md5}
compute-md5 {snakemake.output.ped} {snakemake.output.ped_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
{DEF_HELPER_FUNCS}

sleep 1s  # try to wait for log file flush
compute-md5 {snakemake.log.log} {snakemake.log.log_md5}
"""
)

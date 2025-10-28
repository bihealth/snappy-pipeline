from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
export_config = args["config"]

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

# Extract around BED file, if given.
if [[ -n "{export_config[path_exon_bed]}" ]] && [[ "{export_config[path_exon_bed]}" != "None" ]]; then
    set -e
    bcftools view \
        -R {export_config[path_exon_bed]} \
        {snakemake.input.vcf} \
    | bcftools sort -T $TMPDIR \
    | bcftools norm -d all \
    | bgzip -c \
    > $TMPDIR/tmp.vcf.gz
    tabix -f $TMPDIR/tmp.vcf.gz
else
    set -e
    ln -sr {snakemake.input.vcf} $TMPDIR/tmp.vcf.gz
    ln -sr {snakemake.input.vcf}.tbi $TMPDIR/tmp.vcf.gz.tbi
fi

# Execute VarFish Annotator
varfish-annotator \
    annotate \
    -XX:MaxHeapSize=10g \
    -XX:+UseG1GC \
    \
    --release {export_config[release]} \
    \
    --self-test-chr1-only \
    --ref-path {snakemake.input.reference} \
    --db-path {export_config[path_db]} \
    --refseq-ser-path {export_config[path_refseq_ser]} \
    --ensembl-ser-path {export_config[path_ensembl_ser]} \
    --input-ped {snakemake.input.ped} \
    \
    --input-vcf $TMPDIR/tmp.vcf.gz \
    --output-db-info {snakemake.output.db_infos} \
    --output-gts {snakemake.output.gts}

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

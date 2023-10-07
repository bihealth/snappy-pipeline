from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

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

sniffles \
    --threads {snakemake.threads} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --input {snakemake.input.snf} \
    --tandem-repeats {snakemake.config[step_config][sv_calling_wgs][sniffles2][tandem_repeats]} \
    --vcf $TMPDIR/tmp_raw.vcf \
    --threads {snakemake.threads}

# Remove decoy lines as sniffles writes out garbage IDs for some of them.
grep -v '^hs37d5' $TMPDIR/tmp_raw.vcf > $TMPDIR/tmp_nohs37d5.vcf

for snf in {snakemake.input.snf}; do
    sniffles_name=$(basename $snf .snf)
    library_name=$(echo $sniffles_name | rev | cut -d . -f 1 | rev)
    echo "$sniffles_name $library_name" >>$TMPDIR/samples.txt
done

bcftools reheader \
    --samples $TMPDIR/samples.txt \
    $TMPDIR/tmp_nohs37d5.vcf \
| bgzip -c >{snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

# Compute MD5 sums on output files
compute-md5 {snakemake.output.vcf} {snakemake.output.vcf_md5}
compute-md5 {snakemake.output.vcf_tbi} {snakemake.output.vcf_tbi_md5}

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

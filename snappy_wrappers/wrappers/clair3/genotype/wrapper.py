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

#: The model base path comes from the configuration.
model_base_path = snakemake.config["step_config"]["variant_calling"]["clair3"]["model_base_path"]
#: The sequencing platform comes directly from the args.
platform = {
    "ONT": "ont",
    "PacBio": "hifi",
    "Illumina": "ilmn",
}[snakemake.params["args"]["seq_platform"]]
#: The model name has to be derived from the libraryKit aka caller_model of the args.
model_name = None
for entry in snakemake.config["step_config"]["variant_calling"]["clair3"]["model_map"]:
    if entry.model == snakemake.params["args"]["caller_model"]:
        model_name = entry.name
        break
if not model_name:
    raise ValueError("Could not find model for caller_model == {model_name}")

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
# trap "rm -rf $TMPDIR" EXIT

echo -e "chr2\t178525989\t178807423" > $TMPDIR/regions.bed

# Run actual tools --------------------------------------------------------------------------------

# Run Clair3 and generate gVCF files.
run_clair3.sh \
  --bed_fn=$TMPDIR/regions.bed \
  --gvcf \
  --sample_name={snakemake.wildcards.library_name} \
  --bam_fn={snakemake.input.bam} \
  --ref_fn={snakemake.config[static_data_config][reference][path]} \
  --threads={snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][num_threads]} \
  --platform="ont" \                   ## options: {ont,hifi,ilmn}
  --model_path="{model_base_path}/{model_name}" \
  --output=$TMPDIR/out

# Compute MD5 sums on output files
compute-md5 {snakemake.output.gvcf} {snakemake.output.gvcf_md5}
compute-md5 {snakemake.output.gvcf_tbi} {snakemake.output.gvcf_tbi_md5}

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

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
caller_model = snakemake.params["args"]["caller_model"]

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

# Override the "nproc" command in the called Clair3.
mkdir -p $TMPDIR/bin
cat <<"EOF" >$TMPDIR/bin/nproc
#!/usr/bin/bash
echo {snakemake.threads}
EOF
chmod +x $TMPDIR/bin/nproc
export PATH=$TMPDIR/bin:$PATH


# Run actual tools --------------------------------------------------------------------------------

nproc

# Run Clair3 and generate gVCF files.
run_clair3.sh \
  --gvcf \
  --sample_name={snakemake.wildcards.library_name} \
  --bam_fn={snakemake.input.bam} \
  --ref_fn={snakemake.config[static_data_config][reference][path]} \
  --threads={snakemake.config[step_config][variant_calling][clair3][num_threads]} \
  --platform={platform} \
  --model_path="{model_base_path}/{caller_model}" \
  --output=$TMPDIR/out

# Fix the gVCF file and create index
pbgzip -d -c $TMPDIR/out/merge_output.gvcf.gz \
| sed '/^##FORMAT=<ID=GT/i ##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">' \
| sed '/^##FORMAT=<ID=GT/i ##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description="Minimum DP observed within the GVCF block">' \
| pbgzip -c \
> {snakemake.output.gvcf}

tabix -f {snakemake.output.gvcf}

# Copy out the phased VCF files.
cp -a $TMPDIR/out/merge_output.vcf.gz {snakemake.output.vcf}
cp -a $TMPDIR/out/merge_output.vcf.gz.tbi {snakemake.output.vcf_tbi}

# Compute MD5 sums on output files
compute-md5 {snakemake.output.vcf} {snakemake.output.vcf_md5}
compute-md5 {snakemake.output.vcf_tbi} {snakemake.output.vcf_tbi_md5}
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

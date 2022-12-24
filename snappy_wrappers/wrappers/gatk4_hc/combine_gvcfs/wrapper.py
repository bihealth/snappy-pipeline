from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

# Write files for reproducibility -----------------------------------------------------------------

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

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

gatk \
    CombineGVCFs \
    --java-options '-Xmx6g -Djava.io.tmpdir=$TMPDIR' \
    --tmp-dir $TMPDIR \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --output {snakemake.output.gvcf} \
    $(for annotation in {snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][annotations]}; do \
        echo --annotation $annotation; \
    done) \
    $(for annotation_group in {snakemake.config[step_config][variant_calling][gatk4_hc_gvcf][annotation_groups]}; do \
        echo --annotation-group $annotation_group; \
    done) \
    $(for path in {snakemake.input.gvcf}; do \
        echo --variant $path; \
    done)

# Compute MD5 sums on output files
md5sum {snakemake.output.gvcf} >{snakemake.output.gvcf_md5}
md5sum {snakemake.output.gvcf_tbi} >{snakemake.output.gvcf_tbi_md5}

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
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

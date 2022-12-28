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

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Run actual tools --------------------------------------------------------------------------------

basedir=$(dirname $(dirname {snakemake.output.vcf}))
workdir=$basedir/work
outdir=$basedir/out

# Ensure the working directory is removed, configManta.py will bail out if it already exists
# trap "rm -rf \"$workdir\"" EXIT

configManta.py \
    --exome \
    --referenceFasta {snakemake.config[static_data_config][reference][path]} \
    --runDir $workdir \
    $(echo "{snakemake.input}" | tr ' ' '\n' | grep -v 'bai$' | sed 's/^/--bam /g')

perl -p -i -e 's/isEmail = .*/isEmail = False/g' $workdir/runWorkflow.py

python2 $workdir/runWorkflow.py \
    --jobs 16

cp -ra $workdir/results $outdir
rm -rf $workdir

pushd $outdir
tar czf results.tar.gz results
ln -sr results/variants/diploidSV.vcf.gz $(basename {snakemake.output.vcf})
ln -sr results/variants/diploidSV.vcf.gz.tbi $(basename {snakemake.output.vcf_tbi})
ln -sr results/variants/candidateSV.vcf.gz \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz
ln -sr results/variants/candidateSV.vcf.gz.tbi \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz.tbi

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

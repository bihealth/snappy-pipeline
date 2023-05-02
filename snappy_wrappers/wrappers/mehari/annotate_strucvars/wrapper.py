import os

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Optionally get path to coverage VCF file.
coverage_vcf = " ".join(getattr(snakemake.input, "vcf_cov", []))

# Get shortcut to configuration of varfish_export step
step_name = snakemake.params.args["step_name"]
export_config = snakemake.config["step_config"][step_name]
# Get shortcut to "fix_manta_invs.py" postprocessing script
fix_manta_invs = os.path.join(
    os.path.dirname(__file__),
    "fix_manta_invs.py",
)

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

samples=$(cut -f 2 {snakemake.input.ped} | tr '\n' ',' | sed -e 's/,$//g')

# Fix the Manta inversions
i=0
for vcf in {snakemake.input.vcf}; do
    let "i=$i+1"
    num=$(printf %03d $i)

    python3 {fix_manta_invs} \
        --reference-fasta {snakemake.config[static_data_config][reference][path]} \
        --input-vcf $vcf \
        --output-vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf
    bcftools sort -o $TMPDIR/fixed_bnd_to_inv.$num.vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf

    # Add the missing "GT" tag
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    > $TMPDIR/header.gt.txt

    bcftools annotate \
        -h $TMPDIR/header.gt.txt \
        $TMPDIR/fixed_bnd_to_inv.$num.vcf \
        -O z \
        -o $TMPDIR/final_for_import.$num.vcf.gz
    tabix -s1 -b2 -e2 -f $TMPDIR/final_for_import.$num.vcf.gz
done

# Perform Mehari structural variant annotation.
mehari \
    annotate \
    strucvars \
    --path-db {export_config[path_mehari_db]} \
    --path-input-ped {snakemake.input.ped} \
    $(for p in $TMPDIR/final_for_import.*.vcf.gz; do \
        echo --path-input-vcf $p; \
    done) \
    --path-output-tsv >(gzip -c > {snakemake.output.gts})

# Compute MD5 sums on output files
compute-md5 {snakemake.output.db_infos} {snakemake.output.db_infos_md5}
compute-md5 {snakemake.output.gts} {snakemake.output.gts_md5}
compute-md5 {snakemake.output.feature_effects} {snakemake.output.feature_effects_md5}

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

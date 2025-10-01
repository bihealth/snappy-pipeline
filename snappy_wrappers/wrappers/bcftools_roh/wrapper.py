from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

if path_targets := getattr(snakemake.input, "path_targets", ""):
    path_targets = f"--regions-file {path_targets}"
if path_af_file := getattr(snakemake.input, "path_af_file", ""):
    path_af_file = f"--AF-file {path_af_file}"

ignore_hormef = "--ignore-homref" if args.get("ignore_homref", False) else ""
skip_indels = "--skip-indels" if args.get("skip_indels", False) else ""
rec_rate = f"--rec_rate {args['rec_rate']}" if args.get("rec_rate", 0.0) > 0.0 else ""

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

# Perform ROH Calling -----------------------------------------------------------------------------

out={snakemake.output.txt}
raw_out=${{out%.regions.txt.gz}}.raw.txt.gz

bcftools roh \
    {path_targets} {path_af_file} \
    {ignore_homref} {skip_indels} {rec_rate} \
    --output $raw_out \
    --output-type srz \
    {snakemake.input.vcf}

# Cut out text and BED files.
(
    set +o pipefail
    zcat $raw_out \
    | head -n 3
    zcat $raw_out \
    | tail -n +4 \
    | egrep "^RG|^# RG"
) | bgzip -c \
> {snakemake.output.txt}

compute-md5 {snakemake.output.txt} {snakemake.output.txt_md5}

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

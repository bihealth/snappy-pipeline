from snakemake import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})


shell.executable("/bin/bash")

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
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
# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Compute md5 checksum
md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo "$checksum"
}}

set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

export TMPDIR=$(realpath $TMPDIR)

jar=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)

java -jar $jar -T ReadBackedPhasing \
    -R {snakemake.input.reference} \
    -I {snakemake.input.bam} \
    -L {snakemake.input.vcf} \
    --variant {snakemake.input.vcf} \
    -o {snakemake.output.vcf}

md5 {snakemake.input.vcf}
md5 {snakemake.input.vcf}.tbi
"""
)

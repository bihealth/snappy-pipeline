"""CUBI+Snakemake wrapper code for applying the filter list.
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -euo pipefail

set -x

# "Local" TMPDIR as the scripts try to do "rename" across file sests otherwise
export TMPDIR=$(dirname $(dirname {snakemake.output.vcf}))/tmp
mkdir -p $TMPDIR
trap "rm -rf $TMPDIR" EXIT KILL TERM INT HUP

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

ln -sr {snakemake.input.vcf} {snakemake.output.full}
ln -sr {snakemake.input.vcf}.tbi {snakemake.output.full}.tbi
ln -sr {snakemake.input.vcf}.md5 {snakemake.output.full}.md5
ln -sr {snakemake.input.vcf}.tbi.md5 {snakemake.output.full}.tbi.md5

# Add the "PROTECTED" filter to the vcf header if necessary
filter="FILTER='PASS'"
n=$(bcftools view -h {snakemake.input.vcf} | grep -c '^##FILTER=<ID=PROTECTED,Description="' || true)
if [[ $n -ne 0 ]]
then
    filter="$filter | FILTER~'PROTECTED'"
fi

bcftools view --include "$filter" -O z -o {snakemake.output.vcf} {snakemake.input.vcf}
tabix {snakemake.output.vcf}

tar -zcvf {snakemake.output.log} {snakemake.input.logs}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5 && \
    popd
pushd $(dirname {snakemake.output.log}) && \
    md5sum $(basename {snakemake.output.log}) >$(basename {snakemake.output.log}).md5 && \
    popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

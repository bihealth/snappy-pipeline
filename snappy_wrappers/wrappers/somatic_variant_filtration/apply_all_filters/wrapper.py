"""CUBI+Snakemake wrapper code for applying the filter list."""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
# "Local" TMPDIR as the scripts try to do "rename" across file sests otherwise
export TMPDIR=$(dirname $(dirname {snakemake.output.vcf}))/tmp
mkdir -p $TMPDIR
trap "rm -rf $TMPDIR" EXIT KILL TERM INT HUP

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
"""
)

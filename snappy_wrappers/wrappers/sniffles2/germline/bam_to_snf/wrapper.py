from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT ERR

sniffles \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --input {snakemake.input.bam} \
    --tandem-repeats {snakemake.config[step_config][sv_calling_wgs][sniffles2][tandem_repeats]} \
    --vcf $TMPDIR/tmp.vcf \
    --snf {snakemake.output.snf} \
    --threads {snakemake.threads}

bgzip -c $TMPDIR/tmp.vcf >{snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
md5sum $(basename {snakemake.output.snf}) >$(basename {snakemake.output.snf}).md5
popd
"""
)

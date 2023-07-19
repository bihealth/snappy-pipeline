from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT ERR

sniffles \
    --threads {snakemake.threads} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --input {snakemake.input.snf} \
    --tandem-repeats {snakemake.config[step_config][sv_calling_wgs][sniffles2][tandem_repeats]} \
    --vcf $TMPDIR/tmp_raw.vcf \
    --threads {snakemake.threads}

# Remove decoy lines as sniffles writes out garbage IDs for some of them.
grep -v '^hs37d5' $TMPDIR/tmp_raw.vcf > $TMPDIR/tmp_nohs37d5.vcf

for snf in {snakemake.input.snf}; do
    sniffles_name=$(basename {snakemake.input.snf} .snf)
    library_name=$(echo $sniffles_name | rev | cut -d . -f 1 | rev)
    echo "$sniffles_name $library_name" >>$TMPDIR/samples.txt
done

bcftools reheader \
    --samples $TMPDIR/samples.txt \
    $TMPDIR/tmp_nohs37d5.vcf \
| bgzip -c >{snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
popd
"""
)

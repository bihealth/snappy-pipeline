# -*- coding: utf-8 -*-
# Perform the final step of XHMM: genotyping.

from snakemake.shell import shell

shell(
    r"""
set -x

out_dir=$(dirname {snakemake.output.vcf})
tmp_dir=$out_dir/../tmp

mkdir -p $tmp_dir

cat >$tmp_dir/header.txt <<"EOF"
##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect SV">
EOF

echo -e "1e-08\t6\t70\t-3\t1\t0\t1\t3\t1" >$tmp_dir/params.txt

xhmm --genotype \
    -p $tmp_dir/params.txt \
    -r {snakemake.input.center_zscore} \
    -R {snakemake.input.refilter_original} \
    -g {snakemake.input.discover_xcnv} \
    -F {snakemake.config[static_data_config][reference][path]} \
    -v $tmp_dir/xhmm_genotypes.vcf

bgzip -c $tmp_dir/xhmm_genotypes.vcf \
> $tmp_dir/xhmm_genotypes.vcf.gz
tabix -f $tmp_dir/xhmm_genotypes.vcf.gz

bcftools annotate \
    --header-lines $tmp_dir/header.txt \
    $tmp_dir/xhmm_genotypes.vcf.gz \
| awk -F $'\t' '
    BEGIN {{ OFS = FS; }}
    /^#/ {{ print $0; }}
    /^[^#]/ {{ $8 = $8 ";SVMETHOD=XHMMv0.0.0.2016_01_04.cc14e52"; print $0; }}
    ' \
| bgzip -c \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
popd
"""
)

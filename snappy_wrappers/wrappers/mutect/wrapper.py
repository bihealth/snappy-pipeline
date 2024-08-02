# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect: Snakemake wrapper.py

isort:skip_file
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

intervals = " ".join(snakemake.params.args["intervals"])

shell(
    r"""
set -x

export JAVA_HOME=$(dirname $(which mutect_nonfree))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/out

out_base=$TMPDIR/out/$(basename {snakemake.output.vcf} .vcf.gz)

mutect_nonfree \
    -Xmx2g \
    --analysis_type MuTect \
    --input_file:normal {snakemake.input.normal_bam} \
    --input_file:tumor {snakemake.input.tumor_bam} \
    --dbsnp {snakemake.config[static_data_config][dbsnp][path]} \
    --cosmic {snakemake.config[static_data_config][cosmic][path]} \
    --reference_sequence {snakemake.config[static_data_config][reference][path]} \
    --coverage_file $out_base.full.wig.txt \
    --out $out_base.full.out.txt \
    --vcf $out_base.full.vcf \
    $(if [[ -n "{intervals}" ]]; then
        for itv in "{intervals}"; do \
            echo -n "--intervals $itv "; \
        done; \
    fi)

gzip $out_base.full.out.txt
gzip $out_base.full.wig.txt
bgzip $out_base.full.vcf
tabix -f $out_base.full.vcf.gz

bcftools filter -i "INFO/SOMATIC=1" $out_base.full.vcf.gz -O z -o $out_base.vcf.gz
tabix -f $out_base.vcf.gz

pushd $TMPDIR/out && \
    for f in $(basename {snakemake.output.vcf} .vcf.gz).*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv $out_base.* $(dirname {snakemake.output.vcf})
"""
)

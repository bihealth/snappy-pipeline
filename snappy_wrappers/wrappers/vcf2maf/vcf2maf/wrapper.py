# -*- coding: utf-8 -*-
"""Wrapper for running VCF2MAF incl VEP variant annotation
"""

import pprint

from snakemake.shell import shell

params = snakemake.params.args
step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]
reference = snakemake.config["static_data_config"]["reference"]["path"]


pprint.pprint(params)

shell(
    r"""
set -x
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

CONDAVEP=$(which vep) # Now using bioconda VEP instead of anaconda's (?), it used to be: $(which variant_effect_predictor.pl)
CONDAVEP=$(dirname $CONDAVEP)

zcat -f {snakemake.input.vcf} > $TMPDIR/vcf.vcf

# The cache version must be synchronised with the environment requirements for vep, and of course with the vep_data_path in the config
# WARNING- filter_vcf only works with file ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz (vcf2maf limitation)
vcf2maf.pl --input-vcf $TMPDIR/vcf.vcf \
    --output-maf {snakemake.output.maf} \
    --tmp-dir $TMPDIR \
    --tumor-id {params[tumor_id]} \
    --normal-id {params[normal_id]} \
    --vcf-tumor-id {params[tumor_sample]} \
    --vcf-normal-id {params[normal_sample]} \
    --vep-path $CONDAVEP \
    --vep-data {config[vep_data_path]} \
    --ref-fasta {reference} \
    --ncbi-build {config[ncbi_build]} \
    --cache-version {config[cache_version]} \
    --filter-vcf {config[filter_vcf]} # /fast/groups/cubi/projects/biotools/VEP/static_data/ExAC/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz

pushd $(dirname {snakemake.output.maf})
md5sum $(basename {snakemake.output.maf}) > $(basename {snakemake.output.maf}).md5
popd
"""
)

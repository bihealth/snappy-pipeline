# -*- coding: utf-8 -*-
"""Wrapper for running VCF2MAF incl VEP variant annotation
"""

import pprint
import re
import shutil

from pathlib import Path

from snakemake.shell import shell

params = snakemake.params.args
step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]
reference = snakemake.config["static_data_config"]["reference"]["path"]

# Extract cache version: from user config, otherwise from vep version, otherwise from data path
cache_version = None
if not cache_version:
    if "cache_version" in config and config["cache_version"] != "":
        cache_version = config["cache_version"]
if not cache_version:
    vep = Path(shutil.which("vep")).resolve()
    m = re.compile("^ensembl-vep-([0-9]+)\\..*$").match(vep.parent.name)
    if m:
        cache_version = m.group(1)
if not cache_version:
    data = Path(config["vep_data_path"]) / config["species"]
    pattern = re.compile("^[0-9]+)_" + config["ncbi_build"] + "$")
    for data_path in filter(lambda x: x.is_dir() and pattern.match(x.name), data.iterdir()):
        if cache_version:
            raise ValueError("Multiple valid cache versions, cannot choose one")
        cache_version = pattern.match(data_path).group(1)
if not cache_version:
    raise ValueError("No valid cache version")

filter_vcf = ""
if "filter_vcf" in config:
    filter_vcf = config["filter_vcf"]

vep_custom = " --vep-custom " + config["vep_custom"] if config["vep_custom"] else ""

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

CONDAVEP=$(dirname $(which vep))

vep --help
vep --cache --dir_cache {config[vep_data_path]} --show_cache_info

zcat -f {snakemake.input.vcf} > $TMPDIR/vcf.vcf

# The cache version must be synchronised with the environment requirements for vep, and of course with the vep_data_path in the config
vcf2maf.pl --input-vcf $TMPDIR/vcf.vcf \
    --output-maf {snakemake.output.maf} \
    --tmp-dir $TMPDIR \
    --tumor-id {params[tumor_id]} \
    --normal-id {params[normal_id]} \
    --vcf-tumor-id {params[tumor_sample]} \
    --vcf-normal-id {params[normal_sample]} \
    --species {config[species]} \
    --ncbi-build {config[ncbi_build]} \
    --filter-vcf "{filter_vcf}" \
    --vep-path $CONDAVEP \
    --vep-data {config[vep_data_path]} \
    --ref-fasta {reference} \
    --cache-version {cache_version} \
    {vep_custom}

pushd $(dirname {snakemake.output.maf})
md5sum $(basename {snakemake.output.maf}) > $(basename {snakemake.output.maf}).md5
popd
"""
)

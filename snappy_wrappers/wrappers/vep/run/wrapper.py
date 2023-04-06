# -*- coding: utf-8 -*-
"""Wrapper for running VEP variant annotation
"""

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

# Get shortcuts to static data and step configuration
static_config = snakemake.config["static_data_config"]
config = snakemake.config["step_config"]["somatic_variant_annotation"]["vep"]

shell(
    r"""
set -x

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

vep --offline --cache \
    --dir_cache {config[path_dir_cache]} \
    --species {config[species]} --assembly {config[assembly]} --cache_version {config[cache_version]} \
    $(if [[ -n "{config[transcript_db]}" ]]; then \
        echo "--{config[transcript_db]}"
    fi) \
    $(if [[ "{config[pick]}" = "yes" ]]; then \
        echo "--pick"
    fi) \
    --fasta {static_config[reference][path]} \
    --everything --stats_text --force_overwrite --buffer_size 500 --verbose \
    --input_file {snakemake.input.vcf} \
    --output_file {snakemake.output.vcf} --vcf --compress_output bgzip

tabix {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
f=$(basename {snakemake.output.vcf})
md5sum $f > $f.md5
md5sum $f.tbi > $f.tbi.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

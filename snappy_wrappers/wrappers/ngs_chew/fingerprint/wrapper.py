# -*- coding: utf-8 -*-
"""Wrapper for running ``ngs-chew fingerprint``."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

path_ref = snakemake.config["static_data_config"]["reference"]["path"]
if "hg19" in path_ref or "37" in path_ref:
    genome_release = "GRCh37"
else:
    genome_release = "GRCh38"

shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

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

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/{{out,sorted,sort.tmp}}

if [[ "{library_kit}" == "PacBio HiFi" ]]; then
    preset=map-hifi
elif [[ "{library_kit}" == "PacBio CLR" ]]; then
    preset=map-pb
elif [[ "{library_kit}" == ONT* ]]; then
    preset=map-ont
else
    >&2 echo "Unknown library kit {library_kit}"
    exit 1
fi

ngs-chew fingerprint \
    --min-coverage 5 \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --output-fingerprint {snakemake.ouput.npz} \
    --input-bam {snakemake.input.bam}
    --genome-release {genome_release}

pushd $(dirname {snakemake.output.npz})
md5sum $(basename {snakemake.output.npz}) >$(basename {snakemake.output.npz_md5})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

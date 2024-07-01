# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for GATK UnifiedGenotyper: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# The step key to use.
step_key = snakemake.params.step_key
caller_key = snakemake.params.caller_key

# Build annotations and intervals arguments
arg_annotations = " ".join(
    [
        "--annotation {}".format(anno)
        for anno in snakemake.config["step_config"][step_key][caller_key]["annotations"]
    ]
)
arg_intervals = " ".join(
    ["--intervals {}".format(interval) for interval in snakemake.params["args"]["intervals"]]
)
arg_input_files = " ".join(
    ["--input_file {}".format(fname) for fname in snakemake.input if fname.endswith(".bam")]
)
ped_file = [p for p in snakemake.input if p.endswith(".ped")][0]

allow_seq_dict_incompatibility = snakemake.config["step_config"][step_key][caller_key][
    "allow_seq_dict_incompatibility"
]

shell.executable("/bin/bash")

downsample_to_coverage = snakemake.config["step_config"][step_key][caller_key][
    "downsample_to_coverage"
]

shell(
    r"""
set -x

# TODO: Fix link problems of tabix. What the link?
export JAVA_HOME=$(dirname $(which gatk_nonfree))/..
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

# TODO: add through shell.prefix
export TMPDIR=$HOME/scratch/tmp

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Prepare arguments
if [[ "{allow_seq_dict_incompatibility}" == "True" ]]; then
    arg_seq_dict="-U ALLOW_SEQ_DICT_INCOMPATIBILITY"
else
    arg_seq_dict=
fi

# Call GATK through the bioconda wrapper
gatk_nonfree -Xmx6g -Djava.io.tmpdir=$TMPDIR \
    --analysis_type UnifiedGenotyper \
    $arg_seq_dict \
    -glm BOTH \
    --out $TMPDIR/tmp.vcf.gz \
    --reference_sequence {snakemake.config[static_data_config][reference][path]} \
    --sample_ploidy 2 \
    --dbsnp {snakemake.config[static_data_config][dbsnp][path]} \
    --downsample_to_coverage {downsample_to_coverage} \
    {arg_intervals} \
    {arg_input_files}

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin
# Embed pedigree information --------------------------------------------------
source $(dirname {__file__})/../../tools/add_ped_header.sh
snappy-add_ped_header "{ped_file}" "$TMPDIR/tmp.vcf.gz" /dev/stdout \
> {snakemake.output.vcf}

# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

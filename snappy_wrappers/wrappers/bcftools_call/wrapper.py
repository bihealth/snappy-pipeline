# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools call: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ignore_chroms = snakemake.config["step_config"]["variant_calling"].get("ignore_chroms", [])
if ignore_chroms:
    arg_ignore_chroms = "--ignore-chroms " + " ".join(map(repr, ignore_chroms))
else:
    arg_ignore_chroms = ""

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty,logging disabled" >"{snakemake.log.log}"
    fi
fi

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Create shortcut to reference
export REF={snakemake.config[static_data_config][reference][path]}

# Define function for mpileup, to be used through GNU  parallel
# TODO: Snakemake wrappers use "bash -c" instead of temporary script files, thus we have to generate a script here
cat <<"EOF" >$TMPDIR/samtools_mpileup.sh
#!/bin/bash
set -euo pipefail
IFS=$'\n\t'
set -x

samtools mpileup \
    --BCF \
    --fasta-ref $1 \
    --max-depth {snakemake.config[step_config][variant_calling][bcftools][max_depth]} \
    --max-idepth {snakemake.config[step_config][variant_calling][bcftools][max_indel_depth]} \
    --output-tags DP,AD,ADF,ADR,SP,INFO/AD,INFO/ADF,INFO/ADR \
    --per-sample-mF \
    --redo-BAQ \
    --uncompressed \
    -r $2 \
    $(echo {snakemake.input} | tr ' ' '\n' | grep '\.bam$') \
| bcftools call \
    --multiallelic-caller \
    --format-fields GQ,GP \
    --variants-only
EOF
chmod u+x $TMPDIR/samtools_mpileup.sh


# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

# the variant calling is performed in parallel using GNU parallel
snappy-genome_windows \
    --fai-file $REF.fai \
    --window-size {snakemake.config[step_config][variant_calling][bcftools][window_length]} \
    {arg_ignore_chroms} \
| parallel \
    --keep-order \
    --verbose \
    --max-procs {snakemake.config[step_config][variant_calling][bcftools][num_threads]} \
    "$TMPDIR/samtools_mpileup.sh $REF {{}}" \
| snappy-vcf_first_header \
| bcftools norm \
    --fasta-ref $REF \
    --multiallelics -any \
| snappy-vcf_sort $REF.fai \
| bgzip -c \
> {snakemake.output.vcf}

# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

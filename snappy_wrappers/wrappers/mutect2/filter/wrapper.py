# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
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
export TMPDIR=/fast/users/$USER/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Extract orientation stats for all chunks
mkdir -p $TMPDIR/f1r2
rel_path=$(realpath --relative-to=$TMPDIR/f1r2 {snakemake.input.f1r2})
pushd $TMPDIR/f1r2
tar -zxvf ${{rel_path}}
popd

# Create command line list of all f1r2 stats
chunks=$(ls $TMPDIR/f1r2/*.tar.gz)
cmd=$(echo "$chunks" | tr '\n' ' ' | sed -e "s/ *$//" | sed -e "s/ / -I /g")

# Create orientation model
gatk --java-options '-Xms4000m -Xmx8000m' LearnReadOrientationModel \
    -I $cmd \
    -O $TMPDIR/read-orientation-model.tar.gz

# Workaround problem with bcftools merging inserting missing values (.) in MPOS
zcat {snakemake.input.raw} \
    | awk '{{
        if ($0 ~ /;MPOS=(\.|[0-9-])/) {{
            match($0, /(.+)MPOS=([^;]+)(.+)/, arr);
            gsub(/\./, "-2147483648", arr[2]);
            printf "%sMPOS=%s%s\n", arr[1], arr[2], arr[3];
        }} else {{
            print $0;
        }}
    }}' \
    > $TMPDIR/in.vcf

# Filter calls
gatk --java-options '-Xms4000m -Xmx8000m' FilterMutectCalls \
    --reference {reference} \
    --tumor-segmentation {snakemake.input.segments} \
    --contamination-table {snakemake.input.table} \
    --ob-priors $TMPDIR/read-orientation-model.tar.gz \
    --stats {snakemake.input.stats} \
    --variant $TMPDIR/in.vcf \
    --output {snakemake.output.full}

# Index & move to final dest
tabix -f {snakemake.output.full}

# Keep only PASS variants in main output
bcftools view -i 'FILTER="PASS"' -O z -o {snakemake.output.vcf} {snakemake.output.full}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
fn=$(basename {snakemake.output.vcf})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.tbi})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.full})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.full_tbi})
md5sum $fn > $fn.md5
popd
"""
)

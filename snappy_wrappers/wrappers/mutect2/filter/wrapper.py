# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]

segments = (
    " --tumor-segmentation {} ".format(snakemake.input.segments)
    if "segments" in snakemake.input
    else ""
)
table = (
    " --contamination-table {} ".format(snakemake.input.table) if "table" in snakemake.input else ""
)

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi


# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

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
    > $tmpdir/in.vcf

# Filter calls
gatk --java-options '-Xms4000m -Xmx8000m' FilterMutectCalls \
    --reference {reference} \
    {segments} {table} \
    --ob-priors {snakemake.input.orientation} \
    --stats {snakemake.input.stats} \
    --variant $tmpdir/in.vcf \
    --output $tmpdir/out.vcf

# Extract sample names
grep -E '^##tumor_sample=' $tmpdir/out.vcf | sed -e 's/^##tumor_sample=//' > $tmpdir/tumor.lst

# Extract normal sample(s), if present
if grep -q '^##normal_sample=' "$tmpdir/out.vcf"; then
    grep -E '^##normal_sample=' "$tmpdir/out.vcf" | sed -e 's/^##normal_sample=//' > "$tmpdir/normal.lst"
else
    # No normal sample (tumor-only mode) → create empty file
    > "$tmpdir/normal.lst"
fi


# Validate
num_tumor=$(wc -l < $tmpdir/tumor.lst)
num_normal=$(wc -l < $tmpdir/normal.lst)


if [[ $num_tumor -gt 1 ]]; then
    echo "ERROR: More than one tumor sample found (not supported yet)" >&2
    exit 1

fi

if [[ $num_normal -gt 1 ]]; then
    echo "ERROR: More than one normal sample found (not supported yet)" >&2
    exit 1
fi

# Tumor–Normal case
if [[ $num_normal -eq 1 ]]; then
    cat $tmpdir/normal.lst $tmpdir/tumor.lst > $tmpdir/samples.lst
    bcftools view --samples-file $tmpdir/samples.lst \
        --output-type z \
        --output {snakemake.output.full_vcf} \
        $tmpdir/out.vcf

elif [[ $num_normal -eq 0 && $num_tumor -eq 1 ]]; then
    # Tumor-only case
    bgzip -c $tmpdir/out.vcf > {snakemake.output.full_vcf}
fi



tabix {snakemake.output.full_vcf}
# Keep only PASS variants in main output
bcftools view -i 'FILTER="PASS"' -O z -o {snakemake.output.vcf} {snakemake.output.full_vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
fn=$(basename {snakemake.output.vcf})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.vcf_tbi})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.full_vcf})
md5sum $fn > $fn.md5
fn=$(basename {snakemake.output.full_vcf_tbi})
md5sum $fn > $fn.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

reference = snakemake.input.reference

segments = (
    " --tumor-segmentation {} ".format(snakemake.input.segments)
    if getattr(snakemake.input, "segments", None) is not None
    else ""
)
table = (
    " --contamination-table {} ".format(snakemake.input.table)
    if getattr(snakemake.input, "table", None) is not None
    else ""
)

if java_options := args.get("java_options", ""):
    java_options = f"--java-options '{java_options}'"

extra_arguments = " ".join(args.get("extra_arguments", []))

shell.executable("/bin/bash")

ShellWrapper(snakemake).run(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

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
gatk {java_options} FilterMutectCalls \
    --reference {reference} \
    {segments} {table} \
    --ob-priors {snakemake.input.orientation} \
    --stats {snakemake.input.stats} \
    --variant $TMPDIR/in.vcf \
    --output $TMPDIR/out.vcf \
    {extra_arguments}

# Extract sample names
grep -E '^##tumor_sample=' $TMPDIR/out.vcf | sed -e 's/^##tumor_sample=//' > $TMPDIR/tumor.lst

# Extract normal sample(s), if present
if grep -q '^##normal_sample=' "$TMPDIR/out.vcf"; then
    grep -E '^##normal_sample=' "$TMPDIR/out.vcf" | sed -e 's/^##normal_sample=//' > "$TMPDIR/normal.lst"
else
    # No normal sample (tumor-only mode) → create empty file
    > "$TMPDIR/normal.lst"
fi


# Validate
num_tumor=$(wc -l < $TMPDIR/tumor.lst)
num_normal=$(wc -l < $TMPDIR/normal.lst)


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
    cat $TMPDIR/normal.lst $TMPDIR/tumor.lst > $TMPDIR/samples.lst
    bcftools view --samples-file $TMPDIR/samples.lst \
        --output-type z \
        --output {snakemake.output.full_vcf} \
        $TMPDIR/out.vcf

elif [[ $num_normal -eq 0 && $num_tumor -eq 1 ]]; then
    # Tumor-only case
    bgzip -c $TMPDIR/out.vcf > {snakemake.output.full_vcf}
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


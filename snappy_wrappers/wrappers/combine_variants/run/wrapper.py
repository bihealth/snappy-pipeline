# -*- coding: utf-8 -*-
"""Wrapper for combining germline & somatic variants"""

from math import ceil

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
sample_name = args.get("sample_name", "")

if mem_mb := snakemake.resources.get("mem_gb", None):
    mem_mb *= 1024
else:
    mem_mb = snakemake.resources.get("mem_mb", 2048)

java_mem_mb = ceil(mem_mb * 0.75)

ShellWrapper(snakemake).run(
    r"""
tmp=$(mktemp -d)
gatk=$(find ${{CONDA_PREFIX:-}} -name GenomeAnalysisTK.jar)

if [[ -n "{sample_name}" ]]
then
    germline=$tmp/germline.vcf.gz
    bcftools reheader --samples <(echo "{sample_name}") {snakemake.input.germline_vcf} > $germline
    tabix $germline
else
    germline={snakemake.input.germline_vcf}
fi

temp_somatic=$tmp/temp_somatic.vcf.gz
bcftools view \
    --samples-file <(echo "{args[tumor_library]}") \
    --output-type z --output $temp_somatic --write-index=tbi \
    {snakemake.input.somatic_vcf}
somatic=$tmp/somatic.vcf.gz
bcftools annotate \
    --annotations $temp_somatic \
    --mark-sites +SOMATIC \
    --output-type z --output $somatic --write-index=tbi \
    $temp_somatic

java -Xmx{java_mem_mb}m -jar $gatk -T CombineVariants \
    --assumeIdenticalSamples \
    -R {snakemake.input.reference} \
    --variant $germline \
    --variant $somatic \
    {extra_arguments} \
    -o $tmp/combined.vcf.gz

bcftools sort --max-mem {mem_mb}M \
    --temp-dir $tmp/sort \
    --output-type z --output {snakemake.output.vcf} --write-index=tbi \
    $tmp/combined.vcf.gz
"""
)

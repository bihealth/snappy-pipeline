# -*- coding: utf-8 -*-
"""Wrapper for combining germline & somatic variants"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
sample_name = args.get("sample_name", "")

if mem_mb := snakemake.resources.get("mem_gb", None):
    mem_mb *= 1024
else:
    mem_mb = snakemake.resources.get("mem_mb", 2048)

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

somatic=$tmp/somatic.vcf.gz
bcftools view \
    --samples-file <(echo "{args[tumor_library]}") \
    --output-type z --output $somatic --write-index=tbi \
    {snakemake.input.somatic_vcf}

java -Xmx{mem_mb}m -jar $gatk -T CombineVariants \
    --assumeIdenticalSamples \
    -R {snakemake.input.reference} \
    --variant $germline \
    --variant $somatic \
    -o $tmp/combined.vcf.gz

bcftools sort --max-mem {mem_mb}M \
    --temp-dir $tmp/sort \
    --output-type z --output {snakemake.output.vcf} --write-index=tbi \
    $tmp/combined.vcf.gz
"""
)

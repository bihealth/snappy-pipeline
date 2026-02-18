# -*- coding: utf-8 -*-
"""Wrapper for combining germline & somatic variants"""

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
sample_name = args.get("sample_name", "")

if mem_mb := snakemake.resources.get("mem_gb", None):
    mem_mb *= 1024
else:
    mem_mb = snakemake.resources.get("mem_mb", 2048)

shell(
    r"""
set -x

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

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

tmp=$(mktemp -d)
gatk=$(find $CONDA_PREFIX -name GenomeAnalysisTK.jar)

if [[ -n "{sample_name}" ]]
then
    vcf=$tmp/reheaded.vcf.gz
    bcftools reheader --samples <(echo "{sample_name}") {snakemake.input.germline_vcf} > $vcf
    tabix $vcf
else
    vcf={snakemake.input.germline_vcf}
fi

java -Xmx{mem_mb}m -jar $gatk -T CombineVariants \
    -R {snakemake.input.reference} \
    --variants $vcf \
    --variants {snakemake.input.somatic_vcf} \
    -o $tmp/combined.vcf

bcftools sort --max-mem {mem_mb}M \
    --temp-dir $tmp/sort \
    --output-type z --output {snakemake.output.vcf} --write-index=tbi

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

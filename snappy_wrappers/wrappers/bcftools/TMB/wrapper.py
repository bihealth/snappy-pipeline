# -*- coding: utf-8 -*-
"""Wrapper for calculating tumor mutation burde with bcftools"""

import hashlib
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

args = getattr(snakemake.params, "args", {})
target_regions = args["target_regions"]
has_annotation = args["has_annotation"]


def _file_md5(path: str) -> str:
    with open(path, "rb") as inputf:
        return hashlib.md5(inputf.read()).hexdigest()


bed_md5 = _file_md5(target_regions)
vcf_md5 = _file_md5(str(snakemake.input.vcf))

missense_re = args["missense_re"] if has_annotation else ""

ShellWrapper(snakemake).run(
    r"""
# Ensure locale is set to C, such that the printf %f calls work correctly
export LC_ALL=C

bed_file={target_regions}
bed_file_name=$(basename $bed_file)

name_vcf=$(basename {snakemake.input.vcf})

# Avoids script failing with gzip error status
cmd=zcat
gzip -t $bed_file || cmd=cat

total_exom_length=$($cmd $bed_file | \
    awk '{{dis+=$3-$2}} END {{print dis}}') #TMB_rounded=`printf "%.3f" $TMB`

number_snvs=$(bcftools view -R $bed_file -v snps --threads 2 -H {snakemake.input.vcf}| wc -l)
number_indels=$(bcftools view -R $bed_file -v indels --threads 2 -H {snakemake.input.vcf}| wc -l)
number_variants=$(bcftools view -R $bed_file --threads 2 -H {snakemake.input.vcf}| wc -l)

if [[ -n "{missense_re}" ]]
then
    number_missense_variants=$(bcftools view -R $bed_file --threads 2 -H {snakemake.input.vcf} | grep -E '{missense_re}' | wc -l || true)
else
    number_missense_variants=0
fi

TMB=$(printf "%f" $(echo "1000000*($number_variants/$total_exom_length)" | bc -l))
missense_TMB=$(printf "%f" $(echo "1000000*($number_missense_variants/$total_exom_length)" | bc -l))
if [[ $(echo "{has_annotation}" | tr '[a-z]' '[A-Z]') = "TRUE" ]]
then
    cat << EOF > {snakemake.output.json}
{{
    "Library_name": "{snakemake.wildcards.tumor_library}",
    "VCF_file": "$name_vcf",
    "VCF_md5": "{vcf_md5}",
    "BED_file": "$bed_file_name",
    "BED_md5": "{bed_md5}",
    "TMB": $TMB,
    "missense_TMB": $missense_TMB,
    "Number_variants": $number_variants,
    "Number_snvs": $number_snvs,
    "Number_indels": $number_indels,
    "Number_missense": $number_missense_variants,
    "Total_regions_length": $total_exom_length
}}
EOF
else
    cat << EOF > {snakemake.output.json}
{{
    "Library_name": "{snakemake.wildcards.tumor_library}",
    "VCF_file": "$name_vcf",
    "VCF_md5": "{vcf_md5}",
    "BED_file": "$bed_file_name",
    "BED_md5": "{bed_md5}",
    "TMB": $TMB,
    "Number_variants": $number_variants,
    "Number_snvs": $number_snvs,
    "Number_indels": $number_indels,
    "Total_regions_length": $total_exom_length
}}
EOF
fi

"""
)

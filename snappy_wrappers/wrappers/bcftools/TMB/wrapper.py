# -*- coding: utf-8 -*-
"""Wrapper for calculating tumor mutation burde with bcftools
"""

from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

missense_re = (
    snakemake.params.args["missense_re"]
    if "args" in snakemake.params.keys() and "missense_re" in snakemake.params.args.keys()
    else ""
)

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

bed_file={snakemake.config[step_config][tumor_mutational_burden][target_regions]}
bed_file_name=$(basename $bed_file)
bed_md5=$(md5sum $bed_file | awk '{{print $1}}')

name_vcf=$(basename {snakemake.input.vcf})
vcf_md5=$(md5sum {snakemake.input.vcf} | awk '{{print $1}}')

total_exom_length=$(zcat $bed_file | \
    awk '{{dis+=$3-$2}} END {{print dis}}') #TMB_rounded=`printf "%.3f" $TMB`

number_snvs=$(bcftools view -R $bed_file -v snps --threads 2 -H {snakemake.input.vcf}| wc -l)
number_indels=$(bcftools view -R $bed_file -v indels --threads 2 -H {snakemake.input.vcf}| wc -l)
number_variants=$(bcftools view -R $bed_file --threads 2 -H {snakemake.input.vcf}| wc -l)

if [[ -n "{missense_re}" ]]
then
    number_missense_variants=$(bcftools view -R $bed_file --threads 2 -H {snakemake.input.vcf}| grep -E '{missense_re}' | wc -l)
fi

TMB=`echo "1000000*($number_variants/$total_exom_length)" | bc -l `
missense_TMB=`echo "1000000*($number_missense_variants/$total_exom_length)" | bc -l `
if [[ {snakemake.config[step_config][tumor_mutational_burden][has_annotation]} == "TRUE" ]]
then
    out_file=$(cat << EOF 
    {{
    "Library_name": {snakemake.wildcards.tumor_library},
    "VCF_file": $name_vcf,
    "VCF_md5": $vcf_md5,
    "BED_file": $bed_file_name,
    "BED_md5": $bed_md5,
    "TMB": $TMB,
    "missense_TMB": $missense_TMB,
    "Number_variants": $number_variants,
    "Number_snvs": $number_snvs,
    "Number_indels": $number_indels,
    "Total_regions_length": $total_exom_length
    }}
EOF
    )
    echo $out_file > {snakemake.output.json}
else
    out_file=$(cat << EOF 
    {{
    "Library_name": {snakemake.wildcards.tumor_library},
    "VCF_file": $name_vcf,
    "VCF_md5": $vcf_md5,
    "BED_file": $bed_file_name,
    "BED_md5": $bed_md5,
    "TMB": $TMB,
    "Number_variants": $number_variants,
    "Number_snvs": $number_snvs,
    "Number_indels": $number_indels,
    "Total_regions_length": $total_exom_length
    }}
EOF
    )
    echo $out_file > {snakemake.output.json}
fi

pushd $(dirname {snakemake.output.json})
md5sum $(basename {snakemake.output.json}) > $(basename {snakemake.output.json_md5})
"""
)

# Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
"""
)

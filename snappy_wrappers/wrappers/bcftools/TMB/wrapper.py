# -*- coding: utf-8 -*-
"""Wrapper for calculating tumor mutation burde with bcftools
"""

from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

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

#sample_name={snakemake.input.sample_name} #could query from sth else

total_exom_length=$(zcat {snakemake.config[target_coverage_report][path_target_interval_list_mapping][path]} | \
    awk '{dis+=$3-$2} END {print dis}') #bed file should be in params

number_variants=$(bcftools view -R {snakemake.config[target_coverage_report][path_target_interval_list_mapping][path]} -H {snakemake.input.vcf}| wc -l)
#number_snps=$(bcftools view -R {snakemake.input.bed_file} -v snps -H {snakemake.input.vcf}| wc -l)
#number_indels=$(bcftools view -R {snakemake.input.bed_file} -v indels -H {snakemake.input.vcf}| wc -l)
TMB=`echo "1000000*($number_variants/$total_exom_length)" | bc -l `
TMB_rounded=`printf "%.3f" $TMB`
to_json=$(cat <<EOF
{
    #"Sample": $sample_name,
    "TMB": $TMB_rounded,
    "number_variants": $number_variants,
    "number_snps" : $number_snps,
    "number_indels" : $number_indels,
}
EOF
)
echo $to_json > {snakemake.output.json}
md5sum $(basename {snakemake.output.json}) > $(basename {snakemake.output.json_md5})

"""
)

#Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
"""
)
# -*- coding: utf-8 -*-
"""Wrapper for running Oncotator
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

module purge
module load HTSlib/1.2.1-foss-2015a
module load BCFtools/1.2-foss-2015a
module load Oncotator/v1.8.0.0-foss-2015a-Python-2.7.9

# Shortcut to corpus directory (line length limit...)
corpus={snakemake.config[step_config][somatic_variant_annotation][oncotator][path_corpus]}

# Save original sample names
bcftools view -h {snakemake.input.vcf} | tail -n 1 | cut -f 10- | tr '\t' '\n' \
>{snakemake.output.samples}

# Prepare input VCF file for Oncotator ------------------------------------------------

# Create new samples file with TUMOR/NORMAL
echo -e "TUMOR\nNORMAL" > {snakemake.output.fake_samples}

# Create transmogrified VCF file for the input of Oncotator
bcftools filter \
    -r "1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y" \
    {snakemake.input.vcf} \
| bcftools reheader --samples {snakemake.output.fake_samples} \
> {snakemake.output.vcf_onco_in}

# Call Oncotator with VCF output ------------------------------------------------------

# Perform Oncotator annotation (using fake sample names)
oncotator -v -i VCF -o VCF \
    --db-dir $corpus \
    -c $corpus/override_lists/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
    --log_name $(dirname {snakemake.log})/oncotator.vcf.log \
    {snakemake.output.vcf_onco_in} \
    {snakemake.output.tmp_vcf} \
    {snakemake.params.genome}

# Add back the real sample names
bcftools reheader --samples {snakemake.output.samples} {snakemake.output.tmp_vcf} \
| bgzip -c \
>{snakemake.output.vcf}
tabix {snakemake.output.vcf}

# Compute MD5 sums
pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf_md5}) && \
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi_md5}) && \
    popd

# Call Oncotator with MAF output ------------------------------------------------------

# Perform Oncotator annotation (using fake sample names)
oncotator -v -i VCF -o TCGAMAF \
    --db-dir $corpus \
    -c $corpus/override_lists/tx_exact_uniprot_matches.AKT1_CRLF2_FGFR1.txt \
    --log_name $(dirname {snakemake.log})/oncotator.vcf.log \
    {snakemake.output.vcf_onco_in} \
    {snakemake.output.tmp_maf} \
    {snakemake.params.genome}
bgzip -c {snakemake.output.tmp_maf} >{snakemake.output.maf}

# Compute MD5 sums
pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.maf}) >$(basename {snakemake.output.maf_md5}) && \
    popd
"""
)

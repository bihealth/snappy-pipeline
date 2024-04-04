# -*- coding: utf-8 -*-
"""
Wrapper for preparation the vcf file for somatic neoepitope prediction
"""

import os
from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

exp_format = config["preparation"]["format"]
preparation_vcf = os.path.join(
    os.path.dirname(__file__),
    "comb_rna.py",
)
ensemble_id = (
    "--ignore-ensembl-id-version"
    if config["preparation"]["ignore-ensembl-id-version"]
    else ""
)

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
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
# -----------------------------------------------------------------------------
# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

# Getting region of SNVs only
bcftools filter -i 'TYPE="snp"' {snakemake.input.vcf} | bcftools query -f '%CHROM\t%POS0\t%END\n' > $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.bed

# Getting snvs from RNA sequencing data
bcftools mpileup -Ov \
    --annotate FORMAT/AD,FORMAT/DP \
    -R $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.bed \
    -f {snakemake.config[static_data_config][reference][path]} \
    --per-sample-mF \
    --max-depth 4000 \
    --redo-BAQ -Oz \
    -o $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz \
    {snakemake.input.bam}
tabix -f $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz

#Need to add option for RNA vcf as well
python3 -W ignore {preparation_vcf} \
    --genecode {snakemake.config[static_data_config][features][path]} \
    --input_vcf {snakemake.input.vcf} --format {config[preparation][format]} \
    --rna-vcf $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz \
    --expression_file {snakemake.input.expression} \
    --mode {config[preparation][mode]} \
    -e {config[preparation][expression-column]} \
    -s {snakemake.wildcards.tumor_library} \
    -o /dev/stdout \
    "{ensemble_id}" \
| bgzip -c > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

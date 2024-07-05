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
    "--ignore-ensembl-id-version" if snakemake.params.args["ignore_ensembl_id_version"] else ""
)
max_depth = snakemake.params.args["max_depth"]
format = snakemake.params.args["format"]
mode = snakemake.params.args["mode"]
expression_column = (
    f"--expression-file {snakemake.params.args['expression_column']}" if format == "custom" else ""
)

id_column = f"--expression-file {snakemake.params.args['id_column']}" if format == "custom" else ""

rna_vcf = (
    f"--rna-vcf $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz "
    if format == "snappy_custom"
    else ""
)
expression_file = (
    f"--expression-file {snakemake.input.expression} "
    if format == "snappy_custom"
    else config["preparation"]["expression_file"]
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

if [[ {format}=="snappy_custom" ]]; then
    # Getting region of SNVs only
    bcftools filter -i 'TYPE="snp"' {snakemake.input.vcf} | bcftools query -f '%CHROM\t%POS0\t%END\n' > $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.bed

    # Getting snvs from RNA sequencing data
    bcftools mpileup -Ov \
        --annotate FORMAT/AD,FORMAT/DP \
        -R $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.bed \
        -f {snakemake.config[static_data_config][reference][path]} \
        --per-sample-mF \
        --max-depth {max_depth} \
        --redo-BAQ -Oz \
        -o $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz \
        {snakemake.input.bam}
    tabix -f $TMPDIR/{snakemake.wildcards.tumor_library}.tmp.vcf.gz
fi

#Need to add option for RNA vcf as well
python3 {preparation_vcf} \
    --genecode {snakemake.config[static_data_config][features][path]} \
    --input-vcf {snakemake.input.vcf} --format {format} \
    --mode {mode} \
    {expression_column} \
    {id_column} \
    -s {snakemake.wildcards.tumor_library} \
    -o /dev/stdout \
    {expression_file}\
    {rna_vcf}\
    "{ensemble_id}" \
| bgzip -c > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

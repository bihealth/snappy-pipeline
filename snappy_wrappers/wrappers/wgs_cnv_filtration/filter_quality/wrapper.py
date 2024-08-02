# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper for threshold-based filtration of WGS SV calls."""

# TODO: works for trios, singletons, or if only one parent available but NOT FOR MORE COMPLICATED CASES
# TODO: currently only works for ERDS+SV2

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# shell.executable('/bin/bash') # XXX

base_dir = os.path.dirname(os.path.realpath(__file__))

thresholds = snakemake.config["step_config"]["wgs_cnv_filtration"]["thresholds"][
    snakemake.wildcards.thresholds
]

shell(
    r"""
set -x

# Load library with helper functions.
source {base_dir}/../funcs.sh

# Get name and number of index, father, and mother ------------------------------------------------

index={snakemake.wildcards.index_library}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

# Actual filtration -------------------------------------------------------------------------------

case "{snakemake.wildcards.thresholds}" in
    no_filter)
        # Shortcut when we don't filter.
        cp {snakemake.input.vcf} {snakemake.output.vcf}
        cp {snakemake.input.vcf}.tbi {snakemake.output.vcf}.tbi
    ;;
    *)
        # Build base filter expression for the index.
        case "{thresholds[index_pass]}{thresholds[index_pass_de_novo]}" in
            FalseFalse)
                exp="(1 == 1)"
            ;;
            FalseTrue)
                exp="(FORMAT/DFT[${{index_no}}:*] == \"PASS\")"
            ;;
            TrueFalse)
                exp="(FORMAT/FT[${{index_no}}:*] == \"PASS\")"
            ;;
            TrueTrue)
                exp="((FORMAT/FT[${{index_no}}:*] == \"PASS\")"
                exp+=" && (FORMAT/DFT[${{index_no}}:*] == \"PASS\"))"
            ;;
        esac
        # Extend filter expression if parents are present.
        if [[ "{thresholds[parent_pass]}" == "True" ]] && [[ ! -z "$mother_no" ]]; then
            exp+=" && (FORMAT/FT[${{mother_no}}:*] == \"PASS\")"
        fi
        if [[ "{thresholds[parent_pass]}" == "True" ]] && [[ ! -z "$father_no" ]]; then
            exp+=" && (FORMAT/FT[${{father_no}}:*] == \"PASS\")"
        fi

        # Perform filtration on the FT/DFT using the expression built above.
        bcftools view \
            -i "$exp" \
            -O z \
            -o {snakemake.output.vcf} \
            {snakemake.input.vcf}

        tabix -f {snakemake.output.vcf}
    ;;
esac

# Compute checksums -------------------------------------------------------------------------------

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

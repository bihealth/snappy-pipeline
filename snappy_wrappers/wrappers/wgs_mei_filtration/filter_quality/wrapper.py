# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for quality filter for wgs_mei_filtration.
"""

# TODO: works for trios, singletons, or if only one parent available but NOT FOR MORE COMPLICATED CASES

import os

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
set -x

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={snakemake.wildcards.index_library}
father=$(awk '($2 == "'$index'") {{ print $3; }}' {snakemake.input.ped})
mother=$(awk '($2 == "'$index'") {{ print $4; }}' {snakemake.input.ped})

index_no=$(get_index {snakemake.input.vcf} "$index")
father_no=$(get_index {snakemake.input.vcf} "$father")
mother_no=$(get_index {snakemake.input.vcf} "$mother")

# Perform the actual filtration
lr_var={snakemake.config[step_config][wgs_mei_filtration][thresholds][conservative][lr_var]}
lr_ref={snakemake.config[step_config][wgs_mei_filtration][thresholds][conservative][lr_ref]}

case "{snakemake.wildcards.thresholds}" in
    conservative*)
        # Build base filter expression for conservative case.
        exp="(LR[${{index_no}}] >= $lr_var)"

        if [[ -n "$father_no" ]]; then
            exp+="&& ("
            exp+="(GT[$father_no] == \"alt\" && LR[$father_no] > $lr_var)"
            exp+="|| (GT[$father_no] == \"ref\" && LR[$father_no] < $lr_ref)"
            exp+=")"
        fi

        if [[ -n "$mother_no" ]]; then
            exp+="&& ("
            exp+="(GT[$mother_no] == \"alt\" && LR[$mother_no] > $lr_var)"
            exp+="|| (GT[$mother_no] == \"ref\" && LR[$mother_no] < $lr_ref)"
            exp+=")"
        fi

        bcftools view \
            -i "$exp" \
            -O z \
            -o {snakemake.output.vcf} \
            {snakemake.input.vcf}
    ;;
    *)
        cp {snakemake.input.vcf} {snakemake.output.vcf}
    ;;
esac

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

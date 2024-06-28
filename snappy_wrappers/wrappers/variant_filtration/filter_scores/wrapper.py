# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for scores filter for variant_filtration."""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

# Short-circuit in case of performing no filtration
if snakemake.wildcards.scores == "score_all":
    shell(
        r"""
    # Scores set to "score_all", just link out the data.
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
    """
    )
    sys.exit(0)  # everything went well!


# Actual Filtration -------------------------------------------------------------------------------

# Get shortcut to scores set.
if snakemake.wildcards.scores == "all_scores":
    scores = {"require_coding": False, "require_gerpp_gt2": False, "min_cadd": None}
else:
    scores = snakemake.config["step_config"]["variant_filtration"]["score_thresholds"][
        snakemake.wildcards.scores
    ]

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

# Build filter expression.

filter=""

if [[ "{scores[min_cadd]}" != None ]]; then
    test -z "$filter" || filter+=" && "
    filter+="(CADD_SCORE >= {scores[min_cadd]} || CADD_SCORE == \".\")"
fi

if [[ "{scores[require_gerpp_gt2]}" != False ]]; then
    test -z "$filter" || filter+=" && "
    filter+="(GERPP_GT2 = 1)"
fi

if [[ "{scores[require_coding]}" != False ]]; then
    test -z "$filter" || filter+=" && "
    filter+="(ANN ~ \"|MODERATE\|HIGH|\")"
fi

# Perform filtration (or copy file if all are passing).

if [[ "$filter" == "" ]] || [[ "{snakemake.wildcards.scores}" == all_scores ]]; then
    cp -L {snakemake.input.vcf} {snakemake.output.vcf}
else
    bcftools view \
        -O z \
        -o {snakemake.output.vcf} \
        -i "$filter" \
        {snakemake.input.vcf}
fi

# Build index and MD5 sum files.

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

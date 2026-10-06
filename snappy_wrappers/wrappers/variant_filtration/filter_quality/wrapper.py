# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for quality filter for variant_filtration."""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

args = getattr(snakemake.params, "args", {})

# Short-circuit in case of performing no filtration
if args["filter_mode"] == "no_filter":
    shell(
        r"""
    # Thresholds set to "no_filter", just link out the data.
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
    """
    )
    sys.exit(0)  # everything went well!


# Actual Filtration -------------------------------------------------------------------------------

# Get shortcut to threshold set.
thresholds = args["filter_config"][args["filter_mode"]]

# Get more filter expressions.
if not thresholds.get("include_expressions"):
    include_expressions = ""
else:
    include_expressions = (
        "&& ("
        + " && ".join(map(lambda x: "({})".format(x), thresholds.get("include_expressions")))
        + ")"
    )

shell(
    r"""
set -x

# Load library with helper functions.
source {base_dir}/../../wgs_sv_filtration/funcs.sh

# Get name and number of index, father, and mother.
index={args[index_library]}

index_no=$(get_index {snakemake.input.vcf} "$index")

# Perform the filtration.

if [[ "{args[filter_mode]}" == conservative ]]; then
    include="(GQ[$index_no] >= {thresholds[min_gq]})"
    include+=" && ("
    include+="    ((GT[$index_no] == \"het\") && (DP[$index_no] >= {thresholds[min_dp_het]})) "
    include+="    || ((GT[$index_no] == \"hom\") && (DP[0] >= {thresholds[min_dp_hom]}))"
    include+=")"
    include+="{include_expressions}"
    bcftools view \
        -i "$include" \
        -O z \
        -o {snakemake.output.vcf}  \
        {snakemake.input.vcf}
else
    cp -L {snakemake.input.vcf} {snakemake.output.vcf}
fi

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
"""
)

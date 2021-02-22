# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for inheritance filter for variant_filtration.
"""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

# Actual Filtration -------------------------------------------------------------------------------

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

# The index must not be empty, mother and father can be.
if [[ -z "$index_no" ]]; then
    >&2 echo "Found no index in pedigree?"
    exit 1
fi

if [[ "{snakemake.wildcards.inheritance}" == de_novo ]]; then
    # Perform filtration for de novo variants
    #
    # Build the include filter string.
    include="(INHOUSE_COHORT_AC == 1)"
    include+=" && (DN[$index_no] = \"Y\") && (GT[$index_no] == \"0/1\")"
    # Filter string extension for Jannovar de novo-specific filters.
    include+=" && ("
    include+="  FT[$index_no] !~ \"DeNovoInSibling\""
    include+="  && FT[$index_no] !~ \"DeNovoParentAd2\""
    # Also filter for bad het/hom counts in de novos.
    include+="  && FT[$index_no] !~ \"Min\""
    include+=")"
    test -n "$father_no" && include+="&& (GT[$father_no] == \"0/0\")"
    test -n "$mother_no" && include+="&& (GT[$mother_no] == \"0/0\")"
    # Perform filtration.
    bcftools view \
        -i "$include" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
elif [[ "{snakemake.wildcards.inheritance}" == recessive_hom ]]; then
    # Perform filtration for homozygous variants, heterozygous in both father and mother
    #
    # Build the include filter string.
    include="(INHOUSE_COHORT_AC <= 16)"
    include+=" && (INHERITANCE == \"AR\")"
    include+=" && (INHERITANCE_RECESSIVE_DETAIL == \"AR_HOM_ALT\")"
    include+=" && (GT[$index_no] == \"1/1\")"
    test -n "$father_no" && include+=" && (GT[$father_no] == \"0/1\")"
    test -n "$mother_no" && include+=" && (GT[$mother_no] == \"0/1\")"
    # Perform filtration.
    bcftools view \
        -i "$include" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
## TODO: separate mother and father calls?? - easier to combine later on
## OR will be used in filtering step 6 --> then separate both?
## TODO: this is actually a dominant filtration filter. Should be renamed?
elif [[ "{snakemake.wildcards.inheritance}" == dominant ]]; then
    # Perform filtration for de novo variants
    #
    # Build the inclusion filter string (left and rigt hand side separatdly).
    lhs="(INHOUSE_COHORT_AC <= 8)"
    lhs+=" && (GT[$index_no] == \"0/1\")"
    test -n "$father_no" && lhs+=" && (GT[$father_no] == \"0/0\")"
    test -n "$mother_no" && lhs+=" && (GT[$mother_no] == \"0/1\")"
    rhs="(INHOUSE_COHORT_AC <= 8)"
    rhs+="&& (GT[$index_no] == \"0/1\")"
    test -n "$father_no" && rhs+=" && (GT[$father_no] == \"0/1\")"
    test -n "$mother_no" && rhs+=" && (GT[$mother_no] == \"0/0\")"
    # Perform filtration.
    bcftools view \
        -i "($lhs) || ($rhs)" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
else
    links=1
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.tbi} {snakemake.output.tbi}
    ln -sr {snakemake.input.tbi_md5} {snakemake.output.tbi_md5}
fi

if [[ "${{links-0}}" -ne 1 ]]; then
    tabix -f {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
fi
"""
)

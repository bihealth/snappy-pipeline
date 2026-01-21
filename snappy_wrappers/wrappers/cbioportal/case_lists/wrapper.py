# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out all_cases_with_mutation_data.txt
"""

import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import snakemake

args = getattr(snakemake.params, "args", {})


def write_case_list(case_list_args, outfile):
    """Takes a biomedsheet and writes a case list for all samples with DNA sequencing data"""
    with open(outfile, "w") as f:
        s = "\n".join(
            [
                "cancer_study_identifier: {cancer_study_id}",
                "stable_id: {cancer_study_id}_{stable_id}",
                "case_list_name: {name}",
                "case_list_description: {description}",
                "case_list_category: {category}",
                "case_list_ids: {all_samples}",
                "",
            ]
        ).format(
            cancer_study_id=args["__cancer_study_id"],
            stable_id=case_list_args["stable_id"],
            name=case_list_args["name"],
            description=case_list_args["description"],
            category=case_list_args["category"],
            all_samples="\t".join(case_list_args["samples"]),
        )
        f.write(s)


for case_list, filename in snakemake.output.items():
    assert case_list in args.keys()
    write_case_list(args[case_list], filename)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out all_cases_with_mutation_data.txt
"""
import os
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from snakemake.script import snakemake


def write_case_list(args, outfile):
    """Takes a biomedsheet and writes a case list for all samples with DNA sequencing data"""
    category = os.path.basename(outfile).replace(".txt", "")
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
            cancer_study_id=snakemake.params["args"]["__cancer_study_id"],
            stable_id=args["stable_id"],
            name=args["name"],
            description=args["description"],
            category=args["category"],
            all_samples="\t".join(args["samples"]),
        )
        f.write(s)


for case_list, filename in snakemake.output.items():
    assert case_list in snakemake.params.keys()
    write_case_list(snakemake.params[case_list], filename)

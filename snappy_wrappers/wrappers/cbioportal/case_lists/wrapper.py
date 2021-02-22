# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out all_cases_with_mutation_data.txt
"""

import biomedsheets
from biomedsheets import shortcuts


def write_case_list_sequenced(sheets, config, outfile):
    """Takes a biomedsheet and writes a case list for all samples with DNA sequencing data"""
    samples = []
    for sheet in sheets:
        for p in sheet.bio_entities.values():
            for s in p.bio_samples.values():
                if s.extra_infos["isTumor"]:
                    samples.append(s.name)
    samples = list(set(samples))

    with open(outfile, "w") as f:
        s = "\n".join(
            [
                "cancer_study_identifier: {cancer_study_id}",
                "stable_id: {cancer_study_id}_sequenced",
                "case_list_name: Sequenced tumors",
                "case_list_description: Sequenced tumor samples",
                "case_list_category: all_cases_with_mutation_data",
                "case_list_ids: {all_samples}",
            ]
        ).format(
            cancer_study_id=config["step_config"]["cbioportal_export"]["cancer_study_id"],
            all_samples="\t".join(samples),
        )
        f.write(s)


write_case_list_sequenced(snakemake.params.sheet, snakemake.config, snakemake.output.sequenced)

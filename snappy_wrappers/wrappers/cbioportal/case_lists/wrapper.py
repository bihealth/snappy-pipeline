# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out all_cases_with_mutation_data.txt
and all_cases_with_cna_data.txt
"""


def write_case_list(
    sheets,
    config,
    outfile,
    category="all_cases_with_mutation_data",
    stable_id_suffix="sequenced",
    name="Sequenced tumors",
    description="Sequenced tumor samples",
):
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
                "stable_id: {cancer_study_id}_{stable_id_suffix}",
                "case_list_name: {name}",
                "case_list_description: {description}",
                "case_list_category: {category}",
                "case_list_ids: {all_samples}",
            ]
        ).format(
            cancer_study_id=config["step_config"]["cbioportal_export"]["cancer_study_id"],
            stable_id_suffix=stable_id_suffix,
            name=name,
            description=description,
            category=category,
            all_samples="\t".join(samples),
        )
        f.write(s)


write_case_list(
    snakemake.params.sheet,
    snakemake.config,
    snakemake.output.sequenced,
    category="all_cases_with_mutation_data",
    stable_id_suffix="sequenced",
    name="Sequenced tumors",
    description="Sequenced tumor samples",
)
if snakemake.config["step_config"]["cbioportal_export"]["path_copy_number_step"]:
    write_case_list(
        snakemake.params.sheet,
        snakemake.config,
        snakemake.output.cna,
        category="all_cases_with_cna_data",
        stable_id_suffix="cna",
        name="Tumors with CNA data",
        description="Tumor samples with Copy Number Alteration data",
    )

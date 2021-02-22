# -*- coding: utf-8 -*-
"""cbioportal metadata helper functions"""


def create_study_meta_file(config, outfile):
    """Write study meta file, given study metadata"""
    s = "\n".join(
        [
            "type_of_cancer: {type_of_cancer}",
            "cancer_study_identifier: {cancer_study_id}",
            "name: {study_name}",
            "short_name: {study_name_short}",
            "description: {study_description}",
            "add_global_case_list: true",
        ]
    ).format(
        type_of_cancer=config["step_config"]["cbioportal_export"]["type_of_cancer"],
        cancer_study_id=config["step_config"]["cbioportal_export"]["cancer_study_id"],
        study_description=config["step_config"]["cbioportal_export"]["study_description"],
        study_name=config["step_config"]["cbioportal_export"]["study_name"],
        study_name_short=config["step_config"]["cbioportal_export"]["study_name_short"],
    )
    with open(outfile, "w") as out:
        out.write(s)
        out.write("\n")


create_study_meta_file(snakemake.config, snakemake.output.meta_file)

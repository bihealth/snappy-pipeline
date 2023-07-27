# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out tsv meta_clinical_patient.txt
"""


import csv
import json
import os


def write_clinical_patient_tsv(donors):
    """Takes a biomedsheet and writes a clinical patient tsv for cbioportal, see
    https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md#the-patient-file
    for specification
    """

    # Header lines, first item must start with #
    # attribute Display Names
    NAMES = ["#Patient Identifier", "Dummy value"]
    # attribute Descriptions
    DESC = ["#Patient Identifier", "Dummy value"]
    # attribute Datatype
    DATATYPE = ["#STRING", "STRING"]
    # attribute Priority
    PRIORITY = ["#1", "1"]
    # attribute columns
    COLUMNS = ["PATIENT_ID", "DUMMY"]

    with open(snakemake.output.patient, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # write header
        writer.writerow(NAMES)
        writer.writerow(DESC)
        writer.writerow(DATATYPE)
        writer.writerow(PRIORITY)
        writer.writerow(COLUMNS)

        for donor in donors.keys():
            writer.writerow([donor, "UNKNOWN"])


class SampleInfoTMB:
    step = "tumor_mutational_burden"
    name = "TMB"
    description = "Tumor mutational burden computed on CDS regions"
    datatype = "NUMBER"
    priority = "2"
    column = "TMB"

    def __init__(self, config, **kwargs):
        name_pattern = "bwa." + kwargs["somatic_variant_tool"] + ".tmb.{library}"
        self.tpl = os.path.join(
            config["path"], "output", name_pattern, "out", name_pattern + ".json"
        )

    def get_data(self, lib_by_extraction):
        if "DNA" in lib_by_extraction:
            library = lib_by_extraction["DNA"]
            path = self.tpl.format(library=library)
            try:
                with open(path, "r") as f:
                    result = json.load(f)
                return result["TMB"]
            except Exception as e:
                print(
                    "WARNING- error {} occured when extraction TMB for library {}".format(
                        e, library
                    )
                )
        else:
            print("WARNING- no DNA data")
        return ""


def write_clinical_samples_tsv(donors):
    """Takes a biomedsheet and writes a clinical sample tsv for cbioportal, see
    https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md#the-samples-file
    for specification
    """

    sample_info_getters = []
    config = snakemake.config["step_config"]["cbioportal_export"]
    for step, extra_info in config["sample_info"].items():
        if step == "tumor_mutational_burden":
            sample_info_getters.append(
                SampleInfoTMB(
                    extra_info,
                    somatic_variant_tool=config["somatic_variant_calling_tool"],
                )
            )
        else:
            raise Exception("Unknown sample info request")

    # Header lines, first item must start with #
    # attribute Display Names
    NAMES = ["#Patient Identifier", "Sample Identifier"]
    # attribute Descriptions
    DESC = ["#Patient Identifier", "Sample Identifier"]
    # attribute Datatype
    DATATYPE = ["#STRING", "STRING"]
    # attribute Priority
    PRIORITY = ["#1", "1"]
    # attribute columns
    COLUMNS = ["PATIENT_ID", "SAMPLE_ID"]

    for extra_info in sample_info_getters:
        NAMES += [extra_info.name]
        DESC += [extra_info.description]
        DATATYPE += [extra_info.datatype]
        PRIORITY += [extra_info.priority]
        COLUMNS += [extra_info.column]

    with open(snakemake.output.sample, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # write header
        writer.writerow(NAMES)
        writer.writerow(DESC)
        writer.writerow(DATATYPE)
        writer.writerow(PRIORITY)
        writer.writerow(COLUMNS)

        for donor, v in donors.items():
            for sample, vv in v.items():
                row = [donor, sample]
                for extra_info in sample_info_getters:
                    row.append(extra_info.get_data(vv))
                writer.writerow(row)


write_clinical_patient_tsv(snakemake.params)
write_clinical_samples_tsv(snakemake.params)

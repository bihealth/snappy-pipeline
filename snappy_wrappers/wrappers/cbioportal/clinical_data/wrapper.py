# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing cbioportal patient metadata table from biomedsheets
input. Takes a dict from biomedsheets/snappy_pipeline, writes out tsv meta_clinical_patient.txt
"""

import biomedsheets
from biomedsheets import shortcuts

from pprint import pprint
import csv


def write_clinical_patient_tsv(sheets):
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

    patients = []
    for sheet in sheets:
        patients += [bio_entity.name for bio_entity in sheet.bio_entities.values()]

    with open(snakemake.output.patients_tsv, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # write header
        writer.writerow(NAMES)
        writer.writerow(DESC)
        writer.writerow(DATATYPE)
        writer.writerow(PRIORITY)
        writer.writerow(COLUMNS)

        for p in patients:
            writer.writerow([p, "UNKNOWN"])


def write_clinical_samples_tsv(sheets):
    """Takes a biomedsheet and writes a clinical sample tsv for cbioportal, see
    https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md#the-samples-file
    for specification
    """

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

    with open(snakemake.output.samples_tsv, "w") as tsvfile:
        writer = csv.writer(tsvfile, delimiter="\t")
        # write header
        writer.writerow(NAMES)
        writer.writerow(DESC)
        writer.writerow(DATATYPE)
        writer.writerow(PRIORITY)
        writer.writerow(COLUMNS)

        for sheet in sheets:
            for p in sheet.bio_entities.values():
                for s in p.bio_samples.values():
                    if s.extra_infos["isTumor"]:
                        writer.writerow([p.name, s.name])


write_clinical_patient_tsv(snakemake.params.sheet)
write_clinical_samples_tsv(snakemake.params.sheet)

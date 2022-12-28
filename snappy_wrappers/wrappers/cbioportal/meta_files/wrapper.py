# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for cbioportal metadata preparation
See https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md#clinical-data
"""

from glob import glob
import os

target_dir = "work/upload"
source_dir = os.path.abspath(os.path.dirname(__file__))
study_id = snakemake.config["step_config"]["cbioportal_export"]["cancer_study_id"]

# For each template, replace study_id and write to output in work dir
for tpl in glob(os.path.join(source_dir, "meta*.txt")):
    with open(tpl) as infile, open(os.path.join(target_dir, os.path.basename(tpl)), "w") as outfile:
        for line in infile:
            outfile.write(line.replace("__CANCER_STUDY_IDENTIFIER__", study_id))

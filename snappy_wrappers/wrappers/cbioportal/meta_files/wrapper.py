# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for cbioportal metadata preparation
See https://github.com/cBioPortal/cbioportal/blob/master/docs/File-Formats.md#clinical-data
"""

import os

source_dir = os.path.abspath(os.path.dirname(__file__))

args = getattr(snakemake.params, "args", {})

# For each template, replace study_id and write to output in work dir
for fn in snakemake.output:
    gn = os.path.join(source_dir, os.path.basename(fn))
    assert os.path.exists(gn)
    with open(gn, "rt") as infile, open(fn, "wt") as outfile:
        for line in infile:
            outfile.write(line.format(**args))

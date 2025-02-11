# -*- coding: utf-8 -*-
"""Installation of sequenza non-standard packages"""

import os
import sys

from pathlib import Path

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.



base_dir = os.path.abspath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.utils import install_R_packages

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

dest = os.path.dirname(str(snakemake.output.done))
install_R_packages(dest, os.path.join(os.path.dirname(__file__), "R_environment.json"))
Path(snakemake.output.done).touch()

# -*- coding: utf-8 -*-
"""Installation of sequenza non-standard packages"""

import os
from pathlib import Path
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.utils import install_R_package  # noqa: E402

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

dest = os.path.dirname(str(snakemake.output.done))
for package in snakemake.params.packages:
    install_R_package(dest, package["name"], package["repo"])

Path(str(snakemake.output.done)).touch()

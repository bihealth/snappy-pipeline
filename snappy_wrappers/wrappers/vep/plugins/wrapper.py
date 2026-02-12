# -*- coding: utf-8 -*-
"""Wrapper to import VEP plugins (they must all be in the same directory)"""

import os
import requests

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

for plugin in args:
    fn = os.path.dirname(snakemake.output[0], plugin.name + ".pm")
    if plugin.path:
        os.symlink(os.path.abspath(plugin.path), fn)
    else:
        response = requests.get(plugin.url, stream=True).text
        with open(fn, "wt") as f:
            f.write(response)

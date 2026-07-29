# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code to pull docker container"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

container = getattr(snakemake.params, "container")

ShellWrapper(snakemake).run(
    r"""
apptainer pull --name {snakemake.output.container} {container}
"""
)


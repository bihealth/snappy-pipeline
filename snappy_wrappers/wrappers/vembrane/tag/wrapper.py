# -*- coding: utf-8 -*-
"""Wrapper for running vembrane tag"""

import shlex
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
filter_name = args["filter_name"]
expressions = args["expressions"]
expressions_cmd = " ".join(
    f"--tag {tag}={shlex.quote(str(expr))}" for tag, expr in expressions.items()
)
tag_mode = "--tag-mode fail" if args["tag_mode"] == "include" else "--tag-mode pass"
extra = args["extra_args"]

# Actually run the script.
ShellWrapper(snakemake).run(
    r"""
vembrane tag {tag_mode} {extra} {expressions_cmd} {snakemake.input.vcf} | bgzip -c --threads 4 > {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""
)

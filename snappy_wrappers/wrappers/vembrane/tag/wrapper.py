# -*- coding: utf-8 -*-
"""Wrapper for running vembrane tag or filter."""

import shlex
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
mode = args.get("mode", "tag")
extra = args.get("extra_args", "")

aux = args.get("aux", {})
aux_cmd = " ".join(f"--aux {name}={shlex.quote(str(path))}" for name, path in aux.items())

context = args.get("context", [])
if isinstance(context, str):
    context = [context]
context_cmd = " ".join(f"--context {shlex.quote(str(stmt))}" for stmt in context)

context_files = args.get("context_files", [])
if isinstance(context_files, str):
    context_files = [context_files]
context_file_cmd = " ".join(
    f"--context-file {shlex.quote(str(path))}" for path in context_files
)

ontology = args.get("ontology", "")
ontology_cmd = f"--ontology {shlex.quote(str(ontology))}" if ontology else ""

# Actually run the script.
if mode == "tag":
    expressions = args["expressions"]
    expressions_cmd = " ".join(
        f"--tag {tag}={shlex.quote(str(expr))}" for tag, expr in expressions.items()
    )
    tag_mode = "--tag-mode pass"
    cmd = rf"""
vembrane tag {tag_mode} {extra} {expressions_cmd} {snakemake.input.vcf} | bgzip -c --threads 4 > {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""
else:
    expression = shlex.quote(str(args["expression"]))
    cmd = rf"""
vembrane filter {aux_cmd} {context_cmd} {context_file_cmd} {ontology_cmd} {extra} {expression} {snakemake.input.vcf} | bgzip -c --threads 4 > {snakemake.output.vcf}
tabix {snakemake.output.vcf}
"""

ShellWrapper(snakemake).run(cmd)

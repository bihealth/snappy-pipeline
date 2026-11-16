# -*- coding: utf-8 -*-
"""Wrapper for running Mehari variant annotation (mehari 0.51)

The config section's fields become command line options: ``plugins`` one ``--plugin`` per entry,
the transcript lists comma-separated, booleans as flags, unset fields left out.
"""

from snappy_wrappers.snappy_wrapper import ShellWrapper
from snappy_wrappers import utils

__author__ = "Till Hartmann"
__email__ = "till.hartmann@bih-charite.de"

args = getattr(snakemake.params, "args", {})

num_threads = snakemake.threads

mehari_options = " \\\n        ".join(utils.mehari_options(args.get("config", {})))

# handle input database files
tx_dbs = snakemake.input.get("transcripts", [])
if isinstance(tx_dbs, str):
    tx_dbs = [tx_dbs]
tx_args = " ".join([f"--transcripts {tx}" for tx in tx_dbs])

freq_db = snakemake.input.get("frequencies", "")
freq_arg = f"--frequencies {freq_db}" if freq_db else ""

clinvar_db = snakemake.input.get("clinvar", "")
clinvar_arg = f"--clinvar {clinvar_db}" if clinvar_db else ""

ShellWrapper(snakemake).run(
    r"""
set -x

# using --force here because GATK sometimes produces incorrect VCF headers;
# should probably introduce a GATK cleanup rule instead.
bcftools norm --multiallelics -any {snakemake.input.vcf} --threads {num_threads} --force \
    | mehari annotate seqvars \
        {mehari_options} \
        {tx_args} \
        {freq_arg} \
        {clinvar_arg} \
        {extra_arguments} \
        --reference {snakemake.input.reference} \
        --input - \
        --output {snakemake.output.vcf}

# Index the resulting VCF
tabix --threads {num_threads} {snakemake.output.vcf}
"""
)

# -*- coding: utf-8 -*-
"""Wrapper for running Mehari variant annotation (v0.42.0+)"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Till Hartmann"
__email__ = "till.hartmann@bih-charite.de"

args = getattr(snakemake.params, "args", {})
mehari_config = args.get("config", {})

cli_args = []

# skip keys that are handled via snakemake.input
ignore_keys = {"reference", "transcripts", "frequencies", "clinvar"}
num_threads = snakemake.threads

for key, value in mehari_config.items():
    if key in ignore_keys:
        continue

    # skip None, empty strings, and 'none' for sequence reporting
    if value is None or value == "":
        continue
    if key in ("report_cdna_sequence", "report_protein_sequence") and value == "none":
        continue

    # convert snake_case to kebab-case
    kebab_key = key.replace("_", "-")

    # handle boolean flags vs key-value pairs
    if isinstance(value, bool):
        if value:  # Only append the flag if it's True
            cli_args.append(f"--{kebab_key}")
    elif isinstance(value, list):
        for item in value:
            cli_args.append(f"--{kebab_key} {item}")
    else:
        cli_args.append(f"--{kebab_key} {value}")

mehari_options = " \\\n        ".join(cli_args)

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
bcftools norm --multiallelics -any {snakemake.input.vcf} --threads {num_threads} --force | \
  mehari annotate seqvars \
      {mehari_options} \
      {tx_args} \
      {freq_arg} \
      {clinvar_arg} \
      --reference {snakemake.input.reference} \
      --input - \
      --output {snakemake.output.vcf}

# Index the resulting VCF
tabix --threads {num_threads} {snakemake.output.vcf}
"""
)

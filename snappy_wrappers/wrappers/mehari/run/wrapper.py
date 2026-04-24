# -*- coding: utf-8 -*-
"""Wrapper for running Mehari variant annotation (v0.42.0+)"""

from snakemake.shell import shell

__author__ = "Till Hartmann"
__email__ = "till.hartmann@bih-charite.de"

args = getattr(snakemake.params, "args", {})
mehari_config = args.get("config", {})

cli_args = []

# skip keys that are handled via snakemake.input
ignore_keys = {"reference", "transcripts", "frequencies", "clinvar"}

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

shell(
    r"""
set -x

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Run Mehari annotation
mehari annotate seqvars \
    {mehari_options} \
    {tx_args} \
    {freq_arg} \
    {clinvar_arg} \
    --reference {snakemake.input.reference} \
    --input {snakemake.input.vcf} \
    --output {snakemake.output.vcf}

# Index the resulting VCF
tabix {snakemake.output.vcf}

# Compute MD5 sums for outputs
pushd $(dirname {snakemake.output.vcf})
f=$(basename {snakemake.output.vcf})
md5sum $f > $f.md5
md5sum $f.tbi > $f.tbi.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

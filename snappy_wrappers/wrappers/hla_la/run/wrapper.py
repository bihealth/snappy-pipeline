# -*- coding: utf-8 -*-
"""Wrapper for actually running HLA-LA"""

import csv
import json
import subprocess

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

HLA-LA.pl --maxThreads {snakemake.threads} --sampleID {args[sample_id]} \
    --workingDir $TMPDIR \
    --customGraphDir $(dirname {snakemake.input.path_graph}) \
    --BAM {snakemake.input.bam}

cp $TMPDIR/{args[sample_id]}/hla/R1_bestguess_G.txt {snakemake.output.txt}

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
popd

# To make the pipeline happy.
md5sum {snakemake.log.log} > {snakemake.log.log}.md5
"""
)

alleles = {}
with open(snakemake.output.txt, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        if int(row["perfectG"]) == 0 or float(row["Q1"]) < args["min_score"]:
            continue
        if row["Locus"] not in alleles:
            alleles[row["Locus"]] = set()
        alleles[row["Locus"]] |= set([row["Allele"]])

for k, v in alleles.items():
    alleles[k] = list(v)

with open(snakemake.output.json, "wt") as f:
    json.dump(alleles, f, indent=4)

with open(str(snakemake.output.json) + ".md5", "wt") as f:
    subprocess.run(["md5sum", snakemake.output.json], stdout=f)

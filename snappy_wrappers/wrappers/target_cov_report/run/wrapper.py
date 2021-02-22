# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for STAR: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell.executable("/bin/bash")

# Extract NGS library name.
library_name = snakemake.wildcards.mapper_lib.split(".", 1)[1]
# Get map from library to kit name and paths to BED files.
ngs_library_to_kit = snakemake.params.ngs_library_to_kit

# Get path to target BED file for coverage report, if possible depending on
# enrichment kit of sample.
kit_name = ngs_library_to_kit.get(library_name, ngs_library_to_kit.get("default", ""))

for item in snakemake.config["step_config"]["ngs_mapping"]["target_coverage_report"][
    "path_target_interval_list_mapping"
]:
    if item["name"] == kit_name:
        path_targets_bed = item["path"]
        break
else:
    path_targets_bed = ""

shell(
    r"""
set -x

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Get sorted targets BED file.
zcat --force "{path_targets_bed}" \
| awk -F $'\t' 'BEGIN {{ OFS = FS; }} ($2 < $3) {{ print; }}' \
> $TMPDIR/targets.tmp.bed

bedtools sort \
    -i $TMPDIR/targets.tmp.bed \
    -faidx {snakemake.config[static_data_config][reference][path]}.genome \
| uniq \
> $TMPDIR/targets.bed

bedtools coverage \
    -a $TMPDIR/targets.bed \
    -b {snakemake.input.bam} \
    -g {snakemake.config[static_data_config][reference][path]}.genome \
    -hist \
    -sorted \
| python $(dirname {__file__})/../../../tools/bam_cov_stats.py \
    --bed-path $TMPDIR/targets.bed \
    --min-cov-warning {snakemake.config[step_config][ngs_mapping][target_coverage_report][min_cov_warning]} \
    --min-cov-ok {snakemake.config[step_config][ngs_mapping][target_coverage_report][min_cov_ok]} \
    --max-coverage {snakemake.config[step_config][ngs_mapping][target_coverage_report][max_coverage]} \
    $(if [[ "{snakemake.config[step_config][ngs_mapping][target_coverage_report][detailed_reporting]}" == "True" ]]; then \
        echo --report dec; \
    fi) \
> {snakemake.output.txt}

pushd $(dirname {snakemake.output.txt}) &&
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5 &&
    popd
"""
)

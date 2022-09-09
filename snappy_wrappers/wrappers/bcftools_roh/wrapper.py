# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for bcftools roh
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

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

# Perform ROH Calling -----------------------------------------------------------------------------

out={snakemake.output.txt}
raw_out=${{out%.regions.txt.gz}}.raw.txt.gz

bcftools roh \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][path_targets]}" != "None" ]]; then
        echo --regions-file "{snakemake.config[step_config][roh_calling][bcftools_roh][path_targets]}"
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][path_af_file]}" != "None" ]]; then
        echo --AF-file "{snakemake.config[step_config][roh_calling][bcftools_roh][path_af_file]}"
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][ignore_homref]}" != "False" ]]; then
        echo --ignore-homref
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][skip_indels]}" != "False" ]]; then
        echo --skip-indels
    fi) \
    $(if [[ "{snakemake.config[step_config][roh_calling][bcftools_roh][rec_rate]}" != "None" ]]; then
        echo --rec-rate "{snakemake.config[step_config][roh_calling][bcftools_roh][rec_rate]}"
    fi) \
    --output $raw_out \
    --output-type srz \
    {snakemake.input.vcf}

# Cut out text and BED files.
(
    set +o pipefail
    zcat $raw_out \
    | head -n 3
    zcat $raw_out \
    | tail -n +4 \
    | egrep "^RG|^# RG"
) | bgzip -c \
> {snakemake.output.txt}

pushd $(dirname {snakemake.output.txt})
md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
"""
)

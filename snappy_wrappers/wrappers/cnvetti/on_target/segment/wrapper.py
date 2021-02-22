# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's call step
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

shell(
    r"""
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# Define some global shortcuts
REF={snakemake.config[static_data_config][reference][path]}
REGIONS={snakemake.config[step_config][somatic_targeted_seq_cnv_calling][cnvetti_on_target][path_target_regions]}

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

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Perform the actual processing

cnvetti cmd segment \
    -vvv \
    --segmentation HaarSeg \
    --haar-seg-l-max 5 \
    --haar-seg-fdr 0.001 \
    --input {snakemake.input.bcf} \
    --output {snakemake.output.targets_bcf} \
    --output-segments {snakemake.output.segments_bcf}

# Compute MD5 checksums

pushd $(dirname "{snakemake.output.targets_bcf}")
md5sum $(basename "{snakemake.output.targets_bcf}") >$(basename "{snakemake.output.targets_bcf}").md5
md5sum $(basename "{snakemake.output.targets_bcf}").csi >$(basename "{snakemake.output.targets_bcf}").csi.md5
md5sum $(basename "{snakemake.output.segments_bcf}") >$(basename "{snakemake.output.segments_bcf}").md5
md5sum $(basename "{snakemake.output.segments_bcf}").csi >$(basename "{snakemake.output.segments_bcf}").csi.md5

"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

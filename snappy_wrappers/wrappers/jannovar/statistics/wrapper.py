# -*- coding: utf-8 -*-
"""Wrapper for running Jannovar variant annotation
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

shell.executable("/bin/bash")

shell(
    r"""
set -x

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

# See the following for the memory-related massaging
#
# http://bugs.java.com/view_bug.do?bug_id=8043516

# Call jannovar statistics
MALLOC_ARENA_MAX=4 \
jannovar \
    statistics \
    -XX:MaxHeapSize=4g \
    -XX:CompressedClassSpaceSize=1024m \
    -i {snakemake.input.vcf} \
    -o {snakemake.output.report} \
    -d {snakemake.config[step_config][variant_calling][jannovar_statistics][path_ser]}

pushd $(dirname {snakemake.output.report}) && \
    md5sum $(basename {snakemake.output.report}) > $(basename {snakemake.output.report}).md5
"""
)

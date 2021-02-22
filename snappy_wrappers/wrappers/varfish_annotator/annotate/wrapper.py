# -*- coding: utf-8 -*-
"""Wrapper for running ``varfish-annotator annotate``."""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

# Get shortcuts to static data and step configuration
static_config = snakemake.config["static_data_config"]
export_config = snakemake.config["step_config"]["variant_export"]

this_file = __file__

# Build list of arguments to pass to Jannovar
annotation_args = []

# TODO: care about case of WGS data
# TODO: properly handle release
# TODO: remove case ID parameter from annotator

shell(
    r"""
set -x

# TODO: remove this again, is for fail early
# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} {snakemake.log.wrapper}

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT


if [[ {snakemake.params.args[is_wgs]} == True ]]; then
    set -e
    bcftools view \
        -R {export_config[path_exon_bed]} \
        {snakemake.input.vcf} \
    | bcftools sort -T $TMPDIR \
    | bcftools norm -d all \
    | bgzip -c \
    > $TMPDIR/tmp.vcf.gz
    tabix -f $TMPDIR/tmp.vcf.gz
else
    set -e
    ln -s {snakemake.input.vcf} $TMPDIR/tmp.vcf.gz
    ln -s {snakemake.input.vcf}.tbi $TMPDIR/tmp.vcf.gz.tbi
fi


# See the following for the memory-related massaging
#
# http://bugs.java.com/view_bug.do?bug_id=8043516

out_db_info={snakemake.output.db_infos}
out_gts={snakemake.output.gts}

# Call jannovar statistics
MALLOC_ARENA_MAX=4 \
varfish-annotator \
    annotate \
    -XX:MaxHeapSize=10g \
    -XX:+UseConcMarkSweepGC \
    \
    --release GRCh37 \
    \
    --ref-path {snakemake.config[static_data_config][reference][path]} \
    --db-path {export_config[path_db]} \
    --refseq-ser-path {export_config[path_refseq_ser]} \
    --ensembl-ser-path {export_config[path_ensembl_ser]} \
    \
    --input-vcf $TMPDIR/tmp.vcf.gz \
    --output-db-info ${{out_db_info%.gz}} \
    --output-gts ${{out_gts%.gz}}

gzip ${{out_db_info%.gz}}
gzip ${{out_gts%.gz}}

pushd $(dirname {snakemake.output.gts}) && \
    md5sum $(basename {snakemake.output.gts}) > $(basename {snakemake.output.gts}).md5 && \
    md5sum $(basename {snakemake.output.db_infos}) > $(basename {snakemake.output.db_infos}).md5
"""
)

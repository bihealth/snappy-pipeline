# -*- coding: utf-8 -*-
"""Wrapper for running ``varfish-annotator annotate``."""

import os

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

# Get shortcuts to static data and step configuration
static_config = snakemake.config["static_data_config"]
export_config = snakemake.config["step_config"][snakemake.params.args["step_name"]]

this_file = __file__

fix_manta_invs = os.path.join(
    os.path.dirname(__file__),
    "fix_manta_invs.py",
)

# Build list of arguments to pass to Jannovar
annotation_args = []

# TODO: care about case of WGS data
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
#trap "rm -rf $TMPDIR" EXIT


# See the following for the memory-related massaging
#
# http://bugs.java.com/view_bug.do?bug_id=8043516

out_db_info={snakemake.output.db_infos}
out_gts={snakemake.output.gts}
out_effects={snakemake.output.feature_effects}

samples=$(cut -f 2 {snakemake.input.ped} | tr '\n' ',' | sed -e 's/,$//g')

# Fix the Manta inversions
python3 {fix_manta_invs} \
    --reference-fasta {snakemake.config[static_data_config][reference][path]} \
    --input-vcf {snakemake.input.vcf} \
    --output-vcf $TMPDIR/fixed_bnd_to_inv_unsorted.vcf
bcftools sort -o $TMPDIR/fixed_bnd_to_inv.vcf $TMPDIR/fixed_bnd_to_inv_unsorted.vcf

# TODO: remove vcf4.3 to vcf4.2 conversion once varfish-annotator uses modern HTSJDK
bcftools view \
    --threads 4 \
    --force-samples \
    -s $samples \
    $TMPDIR/fixed_bnd_to_inv.vcf \
| perl -p -i -e 's/fileformat=VCFv4.3/fileformat=VCFv4.2/' \
> $TMPDIR/cut_samples.vcf

echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
> $TMPDIR/header.gt.txt

bcftools annotate \
    -h $TMPDIR/header.gt.txt \
    $TMPDIR/cut_samples.vcf \
| bcftools view \
    -i 'GT ~ "1"' \
    -O z \
    -o $TMPDIR/final_for_import.vcf.gz
tabix -s1 -b2 -e2 -f $TMPDIR/final_for_import.vcf.gz

# Compatibility mode with VarFish Server
compatibility_option=""
if [[ "{snakemake.params.args[varfish_server_compatibility]}" != "False" ]]; then
    compatibility_option="--opt-out callers-array"
fi

# Call jannovar statistics
MALLOC_ARENA_MAX=4 \
varfish-annotator \
    annotate-svs \
    -XX:MaxHeapSize=10g \
    -XX:+UseG1GC \
    \
    --release {export_config[release]} \
    \
    --db-path {export_config[path_db]} \
    --refseq-ser-path {export_config[path_refseq_ser]} \
    --ensembl-ser-path {export_config[path_ensembl_ser]} \
    --input-ped {snakemake.input.ped} \
    \
    --input-vcf $TMPDIR/final_for_import.vcf.gz \
    --output-db-info ${{out_db_info%.gz}} \
    --output-gts ${{out_gts%.gz}} \
    --output-feature-effects ${{out_effects%.gz}}  $compatibility_option

gzip ${{out_db_info%.gz}}
gzip ${{out_gts%.gz}}
gzip ${{out_effects%.gz}}

pushd $(dirname {snakemake.output.gts}) && \
    md5sum $(basename {snakemake.output.feature_effects}) > $(basename {snakemake.output.feature_effects}).md5 && \
    md5sum $(basename {snakemake.output.gts}) > $(basename {snakemake.output.gts}).md5 && \
    md5sum $(basename {snakemake.output.db_infos}) > $(basename {snakemake.output.db_infos}).md5
"""
)

# -*- coding: utf-8 -*-
# Use GATK 3 for computing depth of coverage as input for XHMM.

import os

from snakemake.shell import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.utils import DictQuery

# Pick the target BED file to use.  If it goes by the name "default" then we
# simply take the one at xhmm/path_target_interval_list.  Otherwise, we
# have to go through the list in xhmm/path_target_interval_list_mapping.
xhmm_config = DictQuery(snakemake.config).get("step_config/targeted_seq_cnv_calling/xhmm")
if snakemake.params.args["library_kit"] == "default":
    target_interval_bed = xhmm_config["path_target_interval_list"]
else:
    for item in xhmm_config["path_target_interval_list_mapping"]:
        if item["name"] == snakemake.params.args["library_kit"]:
            target_interval_bed = item["path"]
            break
    else:  # of for, did not break out
        raise Exception("Found no target intervals for %s" % item["name"])

shell(
    r"""
set -x

out_dir=$(dirname {snakemake.output.sample_summary})
out_prefix=$(basename {snakemake.output.sample_summary} .sample_summary)
tmp_dir=$out_dir/../tmp
tmp_list=$tmp_dir/$out_prefix.bam.list

mkdir -p $tmp_dir

# Create BAM list file.
echo $(readlink -f {snakemake.input.bam}) \
> $tmp_list

# Compute coverage.
gatk_nonfree \
    -Xmx3072m \
    -T DepthOfCoverage \
    -I $tmp_list \
    -L {target_interval_bed} \
    -R {snakemake.config[static_data_config][reference][path]} \
    -dt BY_SAMPLE \
    -dcov 5000 \
    -l INFO \
    --omitDepthOutputAtEachBase \
    --omitLocusTable \
    --minBaseQuality 0 \
    --minMappingQuality 20 \
    --start 1 \
    --stop 5000 \
    --nBins 200 \
    --includeRefNSites \
    --countType COUNT_FRAGMENTS \
    -o $out_dir/$out_prefix
"""
)

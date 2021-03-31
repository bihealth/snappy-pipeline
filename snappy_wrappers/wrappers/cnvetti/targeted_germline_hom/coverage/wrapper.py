# -*- coding: utf-8 -*-
# isort:skip_file
import os

from snakemake.shell import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_pipeline.utils import DictQuery

# Pick the target BED file to use.  If it goes by the name "default" then we
# simply take the one at cnvetti_homdel/path_target_interval_list.  Otherwise, we
# have to go through the list in cnvetti_homdel/path_target_interval_list_mapping.
cnvetti_homdel_config = DictQuery(snakemake.config).get(
    "step_config/targeted_seq_cnv_calling/cnvetti_homdel"
)
if snakemake.params.args["library_kit"] == "default":
    target_interval_bed = cnvetti_homdel_config["path_target_interval_list"]
else:
    for item in cnvetti_homdel_config["path_target_interval_list_mapping"]:
        if item["name"] == snakemake.params.args["library_kit"]:
            target_interval_bed = item["path"]
            break
    else:  # of for, did not break out
        raise Exception("Found no target intervals for %s" % item["name"])

shell(
    r"""
set -x

cnvetti cmd coverage \
    --count-kind Fragments \
    --considered-regions TargetRegions \
    --targets-bed {target_interval_bed} \
    --input {snakemake.input.bam} \
    --output {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5
"""
)

# -*- coding: utf-8 -*-
# isort:skip_file

from snakemake.shell import shell

# NOTE: (valid in snappy 0.3 & 0.4 already. Removing the reference to the config in wrappers should not have created the issue.)
#       When called from the sv_calling_wgs, the model has no target interval list,
#       and target_interval_bed is not defined.
#       This is probably incorrect, as intervals should not be pre-processed for WGS data.

args = getattr(snakemake.params, "args", {})
if target_interval_bed := args.get("target_interval_bed", None):
    target_interval_bed = f"-L {target_interval_bed}"

shell(
    r"""
set -x

gatk PreprocessIntervals \
    --bin-length 0 \
    --interval-merging-rule OVERLAPPING_ONLY \
    -R {args[reference]} \
    {target_interval_bed} \
    -O {snakemake.output.interval_list}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.input.reference
args = getattr(snakemake.params, "args", {})

intervals = getattr(snakemake.input, "region", None)
if intervals:
    interval_params = f"--intervals {intervals}"
else:
    interval_params = ""

input_param_list = []
if "normal_bam" in snakemake.input.keys():
    input_param_list.append(f"--input {snakemake.input.normal_bam}")
    if normal := args.get("normal_lib_name", None):
        input_param_list.append(f'--normal "{normal}"')
if "tumor_bam" in snakemake.input.keys():
    input_param_list.append(f"--input {snakemake.input.tumor_bam}")

input_params = " ".join(input_param_list)

if (max_mnp_distance := args.get("max_mnp_distance", -1)) >= 0:
    max_mnp_distance = f"--max-mnp-distance {max_mnp_distance}"
else:
    max_mnp_distance = ""

# The read orientation counts and the stats of the calls, which variant calling declares
if f1r2 := getattr(snakemake.output, "f1r2", None):
    f1r2_param = f"--f1r2-tar-gz {f1r2}"
else:
    f1r2_param = ""
stats = getattr(snakemake.output, "stats", "")

if germline_resource := getattr(snakemake.input, "germline_resource", ""):
    germline_resource_param = f"--germline-resource {germline_resource}"
else:
    germline_resource_param = ""

if panel_of_normals := getattr(snakemake.input, "panel_of_normals", ""):
    panel_of_normals_param = f"--panel-of-normals {panel_of_normals}"
else:
    panel_of_normals_param = ""

if java_options := args.get("java_options", ""):
    java_options = f'--java-options "{java_options}"'

ShellWrapper(snakemake).run(
    r"""
# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

gatk {java_options} Mutect2 \
    --tmp-dir $TMPDIR \
    --reference {reference} \
    {input_params} \
    --output $TMPDIR/calls.vcf \
    {max_mnp_distance} \
    {f1r2_param} \
    {germline_resource_param} \
    {panel_of_normals_param} \
    {interval_params} \
    {extra_arguments}

# Mutect2 writes the stats next to the calls
if [[ -n "{stats}" ]]; then
    mv $TMPDIR/calls.vcf.stats {stats}
fi

bgzip -c $TMPDIR/calls.vcf > {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}
"""
)

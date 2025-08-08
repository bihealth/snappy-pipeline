# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]
config = snakemake.config["step_config"]["somatic_variant_calling"]["mutect2"]
args = snakemake.params.get("args", {})
log = snakemake.log

intervals = getattr(snakemake.input, "region", [])
if intervals:

    def parse_bed_line(line: str) -> tuple[str, int, int]:
        chrom, start, end = line.strip().split("\t")
        return chrom, int(start), int(end)

    with open(intervals, "rt") as f:
        interval_params = " ".join(
            [
                f"--intervals {chrom}:{start + 1}-{end}"
                for chrom, start, end in map(parse_bed_line, f.readlines())
            ]
        )
else:
    interval_params = ""

raw_output = snakemake.output.raw

input_param_list = []
if "normal_bam" in snakemake.input.keys():
    input_param_list.append(
        f'--input {snakemake.input.normal_bam} --normal "{snakemake.params.normal_lib_name}"'
    )
if "tumor_bam" in snakemake.input.keys():
    input_param_list.append(f"--input {snakemake.input.tumor_bam}")

input_params = " ".join(input_param_list)

if (max_mnp_distance := args.get("max_mnp_distance", None)):
    max_mnp_distance = f"--max-mnp-distance {max_mnp_distance}"
else:
    max_mnp_distance = ""

if getattr(snakemake.output, "f1r2", None):
    base_out = str(snakemake.output.raw).removesuffix(".vcf.gz")
    assert str(snakemake.output.f1r2) == base_out + ".f1r2.tar.gz"
    with_f1r2_tar_gz = "yes"
else:
    with_f1r2_tar_gz = ""

if java_options := config.get("java_options"):
    java_options = f'--java-options "{java_options}"'
else:
    java_options = ""

if germline_resource := config.get("germline_resource"):
    germline_resource_param = f"--germline-resource {germline_resource}"
else:
    germline_resource_param = ""

if panel_of_normals := config.get("panel_of_normals"):
    panel_of_normals_param = f"--panel-of-normals {panel_of_normals}"
else:
    panel_of_normals_param = ""

extra_arguments = " ".join(config["extra_arguments"])

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{log.log}" && mkdir -p $(dirname {log.log})
        exec &> >(tee -a "{log.log}" >&2)
    else
        rm -f "{log.log}" && mkdir -p $(dirname {log.log})
        echo "No tty, logging disabled" >"{log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{log.conda_list}
conda info >{log.conda_info}
md5sum {log.conda_list} >{log.conda_list_md5}
md5sum {log.conda_info} >{log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

out_base=$tmpdir/$(basename {raw_output} .vcf.gz)

# Add intervals if required

gatk {java_options} Mutect2 \
    --tmp-dir $tmpdir \
    --reference {reference} \
    {input_params} \
    --output $out_base.vcf \
    {max_mnp_distance} \
    $(if [[ -n "{with_f1r2_tar_gz}" ]]; then \
        echo --f1r2-tar-gz ${{out_base}}.f1r2.tar.gz
    fi) \
    {germline_resource_param} \
    {panel_of_normals_param} \
    {interval_params} \
    {extra_arguments}

rm -f $out_base.vcf.idx

bgzip $out_base.vcf
tabix -f $out_base.vcf.gz

pushd $tmpdir
for f in $out_base.*; do
    md5sum $f >$f.md5
done
popd

mv $out_base.* $(dirname {raw_output})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {log.log} >{log.log_md5}
"""
)

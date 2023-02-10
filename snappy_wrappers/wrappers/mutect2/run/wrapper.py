# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]
config = snakemake.config["step_config"]["somatic_variant_calling"]["mutect2"]

extra_arguments = " ".join(config["extra_arguments"])

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
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

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT
mkdir -p $TMPDIR/out

out_base=$tmpdir/out/$(basename {snakemake.output.raw} .vcf.gz)

# Add intervals if required
intervals=""
if [[ -n "{snakemake.params.args[intervals]}" ]]
then
    for itv in "{snakemake.params.args[intervals]}"
    do
        intervals="$intervals --intervals $itv"
    done
fi

gatk Mutect2 \
    --tmp-dir $tmpdir \
    --input {snakemake.input.normal_bam} \
    --input {snakemake.input.tumor_bam} \
    --normal "{snakemake.params.normal_lib_name}" \
    --reference {reference} \
    --output $out_base.vcf \
    --f1r2-tar-gz ${{out_base}}.f1r2.tar.gz \
    $(if [[ -n "{config[germline_resource]}" ]]; then \
        echo --germline-resource {config[germline_resource]}
    fi) \
    $(if [[ -n "{config[panel_of_normals]}" ]]; then \
        echo --panel-of-normals {config[panel_of_normals]}
    fi) \
    $intervals \
    {extra_arguments}

rm -f $out_base.vcf.idx

bgzip $out_base.vcf
tabix -f $out_base.vcf.gz

# Store the f1r2 tar file in a sub-directory (for compatibility with parallel wrapper)
file_base=$(basename ${{out_base}})
dir_base=$(dirname ${{out_base}})

tar -zcvf ${{out_base}}.f1r2_tar.tar.gz --directory ${{dir_base}} ${{file_base}}.f1r2.tar.gz

pushd $tmpdir
for f in $out_base.*; do
    md5sum $f >$f.md5
done
popd

mv $out_base.* $(dirname {snakemake.output.raw})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)


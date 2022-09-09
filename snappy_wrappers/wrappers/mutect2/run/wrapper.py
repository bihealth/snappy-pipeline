# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MuTect 2: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

reference = snakemake.config["static_data_config"]["reference"]["path"]
germline_resource = snakemake.config["step_config"]["somatic_variant_calling"]["mutect2"][
    "germline_resource"
]
panel_of_normals = snakemake.config["step_config"]["somatic_variant_calling"]["mutect2"][
    "panel_of_normals"
]

shell.executable("/bin/bash")

shell(
    r"""
set -x

# export JAVA_HOME=$(dirname $(which gatk))/..
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# TODO: add through shell.prefix
export TMPDIR=/fast/users/$USER/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/out

out_base=$TMPDIR/out/$(basename {snakemake.output.raw} .vcf.gz)

# Add intervals if required
intervals=""
if [[ -n "{snakemake.params.args[intervals]}" ]]
then
    for itv in "{snakemake.params.args[intervals]}"
    do
        intervals="$intervals --intervals $itv"
    done
fi

# Add panel of normals if required
pon=""
if [[ "X{panel_of_normals}" != "X" ]]
then
    pon=" --panel-of-normals {panel_of_normals}"
fi

gatk Mutect2 \
    --tmp-dir $TMPDIR \
    --input {snakemake.input.normal_bam} \
    --input {snakemake.input.tumor_bam} \
    --normal "{snakemake.params.normal_lib_name}" \
    --reference {reference} \
    --germline-resource {germline_resource} \
    --f1r2-tar-gz ${{out_base}}.f1r2.tar.gz \
    --output $out_base.vcf \
    $intervals $pon

rm -f $out_base.vcf.idx

bgzip $out_base.vcf
tabix -f $out_base.vcf.gz

# Store the f1r2 tar file in a sub-directory (for compatibility with parallel wrapper)
file_base=$(basename ${{out_base}})
dir_base=$(dirname ${{out_base}})

tar -zcvf ${{out_base}}.f1r2_tar.tar.gz --directory ${{dir_base}} ${{file_base}}.f1r2.tar.gz

pushd $TMPDIR && \
    for f in $out_base.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv $out_base.* $(dirname {snakemake.output.raw})
# d=$(dirname {snakemake.output.raw})
# cp /fast/scratch/users/blance_c/tmp/tmpqd1wusadsomatic_variant_calling_mutect2/$d/* $(dirname {snakemake.output.raw})
"""
)

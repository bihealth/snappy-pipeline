# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's call step"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

window_length = args.get("window_length", "")

shell(
    r"""
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# Define some global shortcuts
REF={args[reference]}
REGIONS={args[path_target_regions]}

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
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

# Define bash functions used for processing

cnvetti-coverage()
{{
    set -x
    cnvetti cmd coverage \
        -vvv \
        --reference "$REF" \
        $(if [[ {args[method_name]} == cnvetti_off_target ]]; then
            echo --considered-regions GenomeWide
            echo --mask-piles
            echo --window-length {window_length}
          else
            echo --considered-regions TargetRegions
            echo --targets-bed "$REGIONS"
          fi) \
        --io-threads 1 \
        --count-kind Fragments \
        --input "$1" \
        --output "$2"
}}

cnvetti-normalize()
{{
    set -x
    cnvetti cmd normalize \
        --normalization TotalCoverageSum \
        --input "$1" \
        --output "$2"
}}

cnvetti-merge()
{{
    set -x
    cnvetti cmd merge-cov \
        --output "$3" \
        "$1" \
        "$2"
}}

cnvetti-ratio()
{{
    set -x
    cnvetti cmd ratio \
        --numerator-sample $3 \
        --denominator-sample $4 \
        --output "$2" \
        "$1"
}}

# Perform the actual processing

cnvetti-coverage "{snakemake.input.tumor_bam}" $TMPDIR/cov_tumor.bcf
cnvetti-coverage "{snakemake.input.normal_bam}" $TMPDIR/cov_normal.bcf

cnvetti-normalize $TMPDIR/cov_tumor.bcf $TMPDIR/norm_tumor.bcf
cnvetti-normalize $TMPDIR/cov_normal.bcf $TMPDIR/norm_normal.bcf

cnvetti-merge $TMPDIR/norm_tumor.bcf $TMPDIR/norm_normal.bcf $TMPDIR/norm_merged.bcf

tumor=$(bcftools view -h $TMPDIR/cov_tumor.bcf | grep '^#CHROM' | rev | cut -f 1 | rev)
normal=$(bcftools view -h $TMPDIR/cov_normal.bcf | grep '^#CHROM' | rev | cut -f 1 | rev)

cnvetti-ratio $TMPDIR/norm_merged.bcf "{snakemake.output.bcf}" "$tumor" "$normal"

# Compute MD5 checksums

pushd $(dirname "{snakemake.output.bcf}")
md5sum $(basename "{snakemake.output.bcf}") >$(basename "{snakemake.output.bcf}").md5
md5sum $(basename "{snakemake.output.bcf}").csi >$(basename "{snakemake.output.bcf}").csi.md5
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

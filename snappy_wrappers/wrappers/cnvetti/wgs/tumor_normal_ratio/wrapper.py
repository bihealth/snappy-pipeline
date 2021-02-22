# -*- coding: utf-8 -*-
"""Wrapper for running CNVetti WGS tumor_normal_ratio step.

When a matched normal sample is given for the tumor then a log2-transformed ratio is computed,
otherwise the log2-transformed relative coverage of the tumor is forwarded.
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

tumor_bcf = snakemake.input.tumor_bcf
normal_bcf = getattr(snakemake.input, "normal_bcf", None)

shell(
    r"""
set -x

# Write out information about conda installation --------------------------------------------------

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Setup auto-cleaned TMPDIR -----------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Compute ratio or forward relative coverage ------------------------------------------------------

if [[ "{normal_bcf}" == "None" ]]; then
    cp {tumor_bcf} {snakemake.output.bcf}
    cp {tumor_bcf}.csi {snakemake.output.bcf}.csi
else
    tumor=$(bcftools view -h {tumor_bcf} | grep '^#CHROM' | rev | cut -f 1 | rev)
    normal=$(bcftools view -h {normal_bcf} | grep '^#CHROM' | rev | cut -f 1 | rev)

    cnvetti cmd merge-cov \
        --output $TMPDIR/merged.bcf \
        {tumor_bcf} \
        {normal_bcf}

    cnvetti cmd ratio \
        --numerator-sample $tumor \
        --denominator-sample $normal \
        --output {snakemake.output.bcf} \
        $TMPDIR/merged.bcf

    tabix -f {snakemake.output.bcf}
fi

# Compute MD5 checksums ---------------------------------------------------------------------------

pushd $(dirname "{snakemake.output.bcf}")
md5sum $(basename "{snakemake.output.bcf}") >$(basename "{snakemake.output.bcf}").md5
md5sum $(basename "{snakemake.output.csi}") >$(basename "{snakemake.output.csi}").md5
"""
)

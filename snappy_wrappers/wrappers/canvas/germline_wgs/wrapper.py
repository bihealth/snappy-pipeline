# -*- coding: utf-8 -*-
"""Wrapper for running Canvas in germline variant calling mode on WGS data"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

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

module purge
module load Canvas/1.11.0  # also loads mono
module load HTSlib/1.3.1-foss-2015a

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

mono $EBROOTCANVAS/Canvas.exe Germline-WGS \
    --bam={snakemake.input.bam} \
    --b-allele-vcf={snakemake.input.vcf} \
    --output=$TMPDIR \
    --reference={snakemake.config[step_config][wgs_cnv_calling][canvas][path_reference]} \
    --genome-folder={snakemake.config[step_config][wgs_cnv_calling][canvas][path_genome_folder]} \
    --filter-bed={snakemake.config[step_config][wgs_cnv_calling][canvas][path_filter_bed]} \
    --sample-name={snakemake.wildcards.library_name}

cp $TMPDIR/CNV.vcf.gz {snakemake.output.vcf}

fname={snakemake.output.vcf}
cp $TMPDIR/CNV.CoverageAndVariantFrequency.txt \
    ${{fname%.vcf.gz}}.cov_and_var_freq.txt

tabix -f {snakemake.output.vcf}
pushd $(dirname {snakemake.output.vcf})

for f in *.vcf.gz *.vcf.gz.tbi *.txt; do
    md5sum $f >$f.md5
done
"""
)

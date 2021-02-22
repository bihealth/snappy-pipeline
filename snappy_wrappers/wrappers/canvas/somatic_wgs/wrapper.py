# -*- coding: utf-8 -*-
"""Wrapper for running Canvas in somatic variant calling mode on WGS data
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

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

mono $EBROOTCANVAS/Canvas.exe Somatic-WGS \
    --bam={snakemake.input.tumor_bam} \
    --b-allele-vcf={snakemake.input.somatic_vcf} \
    --output=$(dirname {snakemake.output.vcf}) \
    --reference={snakemake.config[step_config][somatic_wgs_cnv_calling][canvas][path_reference]} \
    --genome-folder={snakemake.config[step_config][somatic_wgs_cnv_calling][canvas][path_genome_folder]} \
    --filter-bed={snakemake.config[step_config][somatic_wgs_cnv_calling][canvas][path_filter_bed]} \
    --sample-name={snakemake.wildcards.cancer_library}

tabix -f {snakemake.output.vcf}
pushd $(dirname {snakemake.output.vcf})

for f in *.vcf.gz *.vcf.gz.tbi; do
    md5sum $f >$f.md5
done
"""
)

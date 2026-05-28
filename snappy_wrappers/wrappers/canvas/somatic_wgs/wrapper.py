# -*- coding: utf-8 -*-
"""Wrapper for running Canvas in somatic variant calling mode on WGS data"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
path_reference = args["path_reference"]
path_genome_folder = args["path_genome_folder"]
path_filter_bed = args["path_filter_bed"]

ShellWrapper(snakemake).run(
    r"""
set -x


module purge
module load Canvas/1.11.0  # also loads mono
module load HTSlib/1.3.1-foss-2015a

mono $EBROOTCANVAS/Canvas.exe Somatic-WGS \
    --bam={snakemake.input.tumor_bam} \
    --b-allele-vcf={snakemake.input.somatic_vcf} \
    --output=$(dirname {snakemake.output.vcf}) \
    --reference={path_reference} \
    --genome-folder={path_genome_folder} \
    --filter-bed={path_filter_bed} \
    --sample-name={args[cancer_library]}

tabix -f {snakemake.output.vcf}
pushd $(dirname {snakemake.output.vcf})

for f in *.vcf.gz *.vcf.gz.tbi; do
    md5sum $f >$f.md5
done
"""
)

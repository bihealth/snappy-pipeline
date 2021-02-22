# -*- coding: utf-8 -*-
"""Wrapper for running SnpEff variant annotation
"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

normal_library = snakemake.params.args["normal_library"]
cancer_library = snakemake.params.args["cancer_library"]

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

module purge
module load Java/1.8.0_92
module load HTSlib/1.3.1-foss-2015a

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

JAR=/fast/users/mholtgr/scratch/build/SnpEff4.3i/snpEff/snpEff.jar

echo -e "{normal_library}\t{cancer_library}" > $TMPDIR/samples.txt

java -Xmx15g -jar $JAR ann \
        -verbose \
        -cancer \
        -cancerSamples $TMPDIR/samples.txt \
        -csvStats {snakemake.output.report_csv} \
        -htmlStats {snakemake.output.report_html} \
        -nodownload \
        -dataDir $(dirname {snakemake.config[step_config][somatic_variant_annotation][snpeff][path_db]}) \
        $(basename {snakemake.config[step_config][somatic_variant_annotation][snpeff][path_db]}) \
        {snakemake.input.vcf} \
> $TMPDIR/$(basename {snakemake.output.vcf} .gz)

bgzip -c $TMPDIR/$(basename {snakemake.output.vcf} .gz) \
> {snakemake.output.vcf}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi}).md5
"""
)

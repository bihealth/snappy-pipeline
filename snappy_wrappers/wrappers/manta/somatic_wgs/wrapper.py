# -*- coding: utf-8 -*-
"""Wrapper for running Manta in somatic variant calling mode on WGS data"""

import os

from snappy_wrappers.simple_wrapper import SimpleWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

outdir = os.path.dirname(snakemake.output.vcf)
workdir = os.path.join(os.path.dirname(outdir), "work")

# Create manta config file if extra_config dict in arguments.
if args.get("extra_config", None) is not None:
    with open(os.path.join(outdir, "config.ini"), "wt") as f:
        print("[manta]", file=f)
        for k, v in args["extra_config"].items():
            print(f"{k} = {v}", file=f)
    config = f"--config {outdir}/config.ini"
else:
    config = ""

cmd=r"""
# configManta.py doesn't like when its output script is already present
rm {workdir}/runWorkflow.py

configManta.py \
    --referenceFasta {snakemake.input.reference} \
    --runDir {workdir} \
    {normal_bam} \
    --tumorBam {snakemake.input.tumor_bam} \
    {exome} {rna} {unstranded} \
    {callRegions} \
    {config}

# Fix issue https://github.com/Illumina/manta/issues/293
sed -ie "s/smtplib.SMTP('localhost')/smtplib.SMTP('localhost', timeout=5)/" \
    {workdir}/runWorkflow.py

python {workdir}/runWorkflow.py \
    -m local \
    -j 16

ln -sr {workdir}/results/variants/somaticSV.vcf.gz {snakemake.output.vcf}
ln -sr {workdir}/results/variants/somaticSV.vcf.gz.tbi {snakemake.output.vcf_tbi}
ln -sr {workdir}/results/variants/candidateSV.vcf.gz {snakemake.output.candidate}
ln -sr {workdir}/results/variants/candidateSV.vcf.gz.tbi {snakemake.output.candidate_tbi}
""".format(
    snakemake=snakemake,
    workdir=workdir,
    normal_bam=f"--normalBam {snakemake.input.normal_bam}" if getattr(snakemake.input, "normal_bam", None) is not None else "",
    exome="--exome" if args.get("exome", False) else "",
    rna="--rna" if args.get("rna", False) else "",
    unstranded="--unstrandedRNA" if args.get("unstranded", False) else "",
    callRegions=f"--callRegions {snakemake.input.callRegions}" if getattr(snakemake.input, "callRegions", None) is not None else "",
    config=config,
)

SimpleWrapper(snakemake).run_bash(cmd)
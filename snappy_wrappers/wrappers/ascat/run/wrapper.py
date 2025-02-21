# -*- coding: utf-8 -*-
"""
Wrapper for actually running ASCAT
Following usage example for version 3.2.0 (commit 3197356)

Mandatory snakemake.input: GCcontent, replictiming, tumor_baf, tumor_logr
Optional snakemake.input: normal_baf, normal_logr

Mandatory snakemake.params.args: chrom_names, gamma, gender, genomeVersion, max_ploidy,
    max_purity, min_ploidy, min_purity, psi_manual, penalty, rho_manual, seed, y_limit
Optional snakemake.params.args: platform (mandatory when no normal)

Mandatory snakemake.output: circos, goodness_of_fit, na, nb, ploidy,
    purity, RData, segments, segments_raw
Optional snakemake.output: 
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.snappy_wrapper import RWrapper

args = getattr(snakemake.params, "args", {})

# TODO: 
cmd=r"""
library(ASCAT)
library(rlang)

plot_dir <- file.path(dirname(dirname("{snakemake.output.segments}")), "plots")
dir.create(plot_dir, recursive=TRUE, mode="0750")

ascat.bc <- ascat.loadData(
    Tumor_LogR_file="{snakemake.input.tumor_logr}",
    Tumor_BAF_file="{snakemake.input.tumor_baf}",
    Germline_LogR_file="{snakemake.input.normal_logr}",
    Germline_BAF_file="{snakemake.input.normal_baf}",
    chrs={chrs},
    gender="{args[gender]}",
    genomeVersion="{args[genomeVersion]}",
    isTargetedSeq=FALSE
)

ascat.plotRawData(ascat.bc, img.dir=plot_dir, img.prefix="Before_correction_")

ascat.bc <- ascat.correctLogR(
    ascat.bc,
    GCcontentfile="{snakemake.input.GCcontent}",
    replictimingfile="{snakemake.input.reptiming}"
)

ascat.plotRawData(ascat.bc, img.dir=plot_dir, img.prefix="After_correction_")

ascat.bc <- ascat.aspcf(
    ASCATobj=ascat.bc,
    penalty={args[penalty]},
    out.dir=dirname("{snakemake.output.RData}"),
    out.prefix="",
    seed={args[seed]}
)

ascat.plotSegmentedData(ascat.bc, img.dir=plot_dir, img.prefix="Segmented_")

ascat.output <- ascat.runAscat(
    ASCATobj=ascat.bc,
    gamma={args[gamma]},
    y_limit={args[y_limit]},
    circos=file.path(plot_dir, "circos"),
    min_ploidy={args[min_ploidy]}, max_ploidy={args[max_ploidy]},
    min_purity={args[min_purity]}, max_purity={args[max_purity]},
    rho_manual={args[rho_manual]}, psi_manual={args[psi_manual]},
    img.dir=plot_dir, img.prefix="ASCAT_",
    write_segments=TRUE
)

QC <- ascat.metrics(ascat.bc, ascat.output)

save(ascat.bc, ascat.output, QC, file="{snakemake.output.RData}")

# Write out results.
write.table(ascat.output$nA, "{snakemake.output.na}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$nB, "{snakemake.output.nb}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$goodnessOfFit, "{snakemake.output.goodness_of_fit}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$ploidy, "{snakemake.output.ploidy}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$segments, "{snakemake.output.segments}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$segments_raw, "{snakemake.output.segments_raw}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$purity, "{snakemake.output.purity}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
""".format(
    snakemake=snakemake,
    args=args,
    chrs='c("{}")'.format('", "'.join(args["chrom_names"])),
)

RWrapper(snakemake).run(cmd)
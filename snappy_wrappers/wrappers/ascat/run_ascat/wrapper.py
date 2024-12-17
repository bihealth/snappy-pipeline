# -*- coding: utf-8 -*-
"""
Wrapper for actually running ASCAT

Mandatory snakemake.input: GCcontent, replictiming, tumor_baf, tumor_logr
Optional snakemake.input: normal_baf, normal_logr

Mandatory snakemake.params.args: gamma, gender, genomeVersion, max_ploidy, max_purity,
    min_ploidy, min_purity, psi_manual, penalty, rho_manual, seed, y_limit
Optional snakemake.params.args: platform (mandatory when no normal)

Mandatory snakemake.output: aberrantcellfraction, circos, goodness_of_fit, na, nb,
    ploidy, purity, RData, segments, segments_raw
Optional snakemake.output: 
"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

from snappy_wrappers.snappy_wrapper import RWrapper

args = getattr(snakemake.params, "args", {})

has_normal = getattr(snakemake.input, "germline_logr", None) is not None and getattr(snakemake.input, "germline_baf", None) is not None
if has_normal:
    germline_genotypes = "NULL"
else:
    germline_genotypes = 'ascat.predictGermlineGenotypes(ascat.bc, platform="{}", img.dir=".", img.prefix=".")'.format(args["platform"])

cmd=r"""
library(ASCAT)

ascat.bc <- ascat.loadData(
    Tumor_LogR_file="{snakemake.input.tumor_logr}",
    Tumor_BAF_file="{snakemake.input.tumor_baf}",
    Germline_LogR_file={normal_logr},
    Germline_BAF_file={normal_baf},
    chrs={chrom_names},
    gender=args[gender],
    genomeVersion=args[genomeVersion],
    isTargetedSeq=FALSE
)

ascat.plotRawData(ascat.bc, img.dir="work", img.prefix="Before_correction_")

ascat.bc <- ascat.correctLogR(
    ascat.bc,
    GCcontentfile="{snakemake.input.GCcontent}",
    replictimingfile="{snakemake.input.reptiming}"
)

ascat.plotRawData(ascat.bc, img.dir="work", img.prefix="After_correction_")

ascat.bc <- ascat.aspcf(
    ASCATobj=ascat.bc,
    ascat.gg={germline_genotypes},
    penalty=args[penalty],
    out.dir=".",
    out.prefix="",
    seed=args[seed]
)

ascat.plotSegmentedData(ascat.bc, img.dir="work", img.prefix="")

ascat.output <- ascat.runAscat(
    ASCATobj=ascat.bc,
    gamma=args[gamma],
    y_limit=args[y_limit],
    circos="{snakemake.output.circos}",
    min_ploidy=args[min_ploidy], max_ploidy=args[max_ploidy],
    min_purity=args[min_purity], max_purity=args[max_purity],
    rho_manual=args[rho_manual], psi_manual=args[psi_manual],
    img.dir=".", img.prefix="",
    write_segments=TRUE
)

QC <- ascat.metrics(ascat.bc, ascat.output)

save(ascat.bc, ascat.output, QC, file="{snakemake.output.RData}")

# Write out results.
write.table(ascat.output$nA, "{snakemake.output.na}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$nB, "{snakemake.output.nb}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$goodnessOfFit, "{snakemake.output.goodness_of_fit}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$ploidy, "{snakemake.output.ploidy}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$psi, "{snakemake.output.psi}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$segments, "{snakemake.output.segments}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$segments_raw, "{snakemake.output.segments_raw}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
write.table(ascat.output$purity, "{snakemake.output.purity}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
""".format(
    snakemake=snakemake,
    args=args,
    germline_logr='"{}"'.format(snakemake.input.germline_logr) if has_normal else "NULL",
    germline_baf='"{}"'.format(snakemake.input.germline_baf) if has_normal else "NULL",
    chrs='c("{}")'.format('", "'.join(args["chrom_names"])),
    germline_genotypes=germline_genotypes,
)

RWrapper(snakemake).run(cmd)
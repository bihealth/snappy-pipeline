# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

cat << _EOF | R --vanilla --slave
segments <- read.table("{snakemake.input.segment}", sep="\t", header=1, stringsAsFactors=FALSE)
stopifnot(all(c("chromosome", "start", "end", "log2", "probes") %in% colnames(segments)))
calls <- read.table("{snakemake.input.call}", sep="\t", header=1, stringsAsFactors=FALSE)
stopifnot(all(c("chromosome", "start", "end", "cn") %in% colnames(calls)))

# Trim segments to avoid edge effects with calls
segments[,"start"] <- segments[["start"]]+1
segments[,"end"]   <- segments[["end"]]-1
segments <- segments[segments[["end"]]-segments[["start"]]>0,]

r1 <- GenomicRanges::GRanges(
    seqnames=segments[["chromosome"]],
    ranges=IRanges::IRanges(start=segments[["start"]], end=segments[["end"]]),
    strand="*"
)
r2 <- GenomicRanges::GRanges(
    seqnames=calls[["chromosome"]],
    ranges=IRanges::IRanges(start=calls[["start"]], end=calls[["end"]]),
    strand="*"
)

i <- GenomicRanges::findOverlaps(r1, r2, ignore.strand=TRUE)
stopifnot(!any(duplicated(S4Vectors::queryHits(i))))

# Default value: diploid segment (number of copies = 2)
segments[,"cn"] <- 2
segments[S4Vectors::queryHits(i),"cn"] <- round(calls[S4Vectors::subjectHits(i),"cn"])

# Reset segment boundaries to original values
segments[,"start"] <- segments[["start"]]-1
segments[,"end"]   <- segments[["end"]]+1

# Add library name
segments[,"ID"] <- "{snakemake.wildcards[library_name]}"

# Rename columns to follow DNAcopy format (as implemented in PureCN)
dna_copy_columns <- c(ID="ID", chromosome="chrom", start="loc.start", end="loc.end", probes="num.mark", log2="seg.mean", cn="C")
segments <- segments[,colnames(segments) %in% names(dna_copy_columns)]
colnames(segments) <- dna_copy_columns[colnames(segments)]
segments <- segments[,dna_copy_columns]

write.table(segments, file="{snakemake.output.final}", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
_EOF
"""
)

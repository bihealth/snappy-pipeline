#' Compute discrete copy number & log2 on genes
#'
#' @param fn names of the cns file
#' @param tx_obj Bioconductor feature description object or its name (human by default, either TxDb or EnsDb), or GTF/GFF features file
#'
#' @example
cns_to_cna <- function(fn, tx_obj=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, pipeline_id="ENSEMBL") {
    if (missing(fn) || is.null(fn) || !is.character(fn) || length(fn)!=1 || is.na(fn) || fn=="")
        stop("Missing or illegal cns filename")

    tx_obj <- get_tx_object(tx_obj)
    genes <- GenomicFeatures::genes(tx_obj)
    stopifnot("gene_id" %in% colnames(GenomicRanges::mcols(genes)))

    segments <- read.table(fn, sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
    stopifnot(all(c("chrom", "loc.start", "loc.end", "seg.mean", "C") %in% colnames(segments)))

    seg <- GenomicRanges::makeGRangesFromDataFrame(segments, seqnames.field="chrom", start.field="loc.start", end.field="loc.end", strand="*", keep.extra.columns=FALSE)

    i <- GenomicRanges::findOverlaps(genes, seg)
    i <- split(S4Vectors::subjectHits(i), S4Vectors::queryHits(i))
    i <- unlist(i[lengths(i) == 1])

    tbl <- data.frame(FeatureID=as.character(genes$gene_id), cn=NA, log2=NA)
    colnames(tbl)[1] <- pipeline_id
    tbl[as.numeric(names(i)),c("cn", "log2")] <- segments[i,c("C", "seg.mean")]

    tbl
}

#' Collapse bins into segments for Control-FREEC output
#'
#' @param ratios Control-FREEC ratio output
#' @param collapse_col column name for segmentation (MedianRatio or CopyNumber)
#'
#' @return a data frame with the chromosome, segment start & end, the collapsed column value & the number of bins in the segment
#'
#' @examples
collapse_ratio <- function(ratios, collapse_col=c("MedianRatio", "CopyNumber"),
                           chrlen=GenomeInfoDb::seqlengths(BSgenome::seqinfo(BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5))) {
    collapse_col <- match.arg(collapse_col)
    # Test chromosome names
    assertthat::assert_that(all(ratios$Chromosome %in% names(chrlen)))
    # On each chromosome...
    y <- lapply(split(ratios, ratios$Chromosome),
                function(x) {
                    # Order rows according to start position
                    x <- x[order(x$Start),]
                    # Add end-of-segment
                    x$End <- c(x$Start[-1]-1, chrlen[x$Chromosome[1]])
                    # Find break points by looking at which differences between successive values are different from 0
                    i <- which(diff(x[,collapse_col])!=0)
                    if (length(i) == 0) {
                        # Single segment
                        return(matrix(c(x$Start[1],
                                        x$End[nrow(x)],
                                        median(x[,collapse_col], na.rm=TRUE),
                                        nrow(x)),
                                      nrow=1,
                                      dimnames=list(NULL, c("start", "end", collapse_col, "bins"))))
                    }
                    # Build an array of row indices of start and end of segment
                    i <- cbind(Start=c(1, i+1), End=c(i, nrow(x)))
                    # Build segment matrix: start & end of segment, segment value and number of bins in segment
                    y <- t(apply(i, 1, function(j) c(x$Start[j[1]],
                                                     x$End[j[2]],
                                                     median(x[j[1]:j[2],collapse_col], na.rm=TRUE),
                                                     j[2]-j[1]+1)))
                    colnames(y) <-  c("start", "end", collapse_col, "bins")
                    y
                })
    do.call("rbind", lapply(names(y), function(i) data.frame(chrom=rep(i, nrow(y[[i]])), y[[i]], stringsAsFactors=F)))
}

#' Convert Control-FREEC ratio data to format suitable for cBioPortal & cnvkit plots
#'
#' The function reads Control-FREEC main file .ratio.txt, and outputs a list with 5 elements
#' the (GISTIC-like) gene calls, the copy-number gene-base value on log2 scale (ratio to ploidity),
#' the segments (chromosome, start, end, mean log2 value (ratio to ploidity) and number of bins), and
#' the similar values in format suitable for cnvkit.
#'
#' Note that the CNV values are on log2 scale. 0 means that a gene or a segment is diploid in the
#' tumor sample. If the tumor sample was entirely made of tumor cells (purity = 1), a CNV value of 1
#' means that there are 4 copies of a gene or a segment in the tumor sample; a CNV value of -1 means
#' that there is only one copy of a gene of a segment in the tumor sample (log2(1/2) = -1).
#' This assumes that the the organism is diploid.
#'
#' @param ratios data.frame containing Control-FREEC ratio data
#' @param ploidity ploidity (default 2)
#' @param small_value offset to avoid infinities on log scale
#' @param org_obj Bioconductor object defining ids mappings in an organism (default human)
#' @param tx_obj Bioconductor TxDb object defining gene locations (default human)
#' @param bs_obj Bioconductor BSgenome object defining the genome sequence(default human)
#' @param GISTIC GISTIC codes
#'
#' @return list with 5 components:
#' Calls- data frame with Hugo_Symbol,
#' Entrez_Gene_Id and GISTIC gene call, Log2- data frame with Hugo_Symbol, Entrez_Gene_id and log2 ratio,
#' Segment- data.frame with chromosome, start & end position, log2 ratio & number of bins
#' CNVkit_cnv- like Segment, with a dummy gene column, and no number of bins
#' CNVkit_cnr- similar to the initial file, only with bins with a median value different from -1
#' @export
#'
#' @examples
freec_to_calls_log2_segments <- function(ratios, ploidity=2, small_value=0.000000001,
                                         org_obj=org.Hs.eg.db::org.Hs.eg.db,
                                         tx_obj=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                         bs_obj=BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5,
                                         GISTIC=data.frame(
                                             Type=c("Deep", "Shallow", "Diploid", "Gain", "Amp"),
                                             GISTIC_value=c(-2, -1, 0, 1, 2),
                                             From=c(0, 1, 2, 3, 6),
                                             To=c(0, 1, 2, 5, Inf),
                                             stringsAsFactors=FALSE
                                         )) {
    if (missing(ratios) || is.null(ratios) || !is.data.frame(ratios) || nrow(ratios)==0 ||
        !all(c("Chromosome", "Start", "CopyNumber", "MedianRatio") %in% colnames(ratios)))
        stop("Missing/illegal Control-FREEC ratio table")
    i <- !is.na(ratios$Chromosome) & as.character(ratios$Chromosome)!="" &
        !is.na(ratios$Start) & !is.na(as.integer(ratios$Start)) & as.integer(ratios$Start)>=0 &
        !is.na(ratios$CopyNumber) & !is.na(as.numeric(ratios$CopyNumber)) &
        !is.na(ratios$CopyNumber) & !is.na(as.numeric(ratios$MedianRatio))
    if (!any(i))
        stop("No valid bin in Control-FREEC ratio table")
    ratios <- ratios[i,,drop=FALSE]
    if (is.null(org_obj) || !is(org_obj, "OrgDb") ||
        !all(c("ENTREZID", "SYMBOL") %in% AnnotationDbi::columns(org_obj)))
        stop("Missing/invalid organism object (Bioconductor OrgDb)")
    if (is.null(tx_obj) || !is(tx_obj, "TxDb"))
        stop("Missing/invalid genomic features object (Bioconductor TxDb)")
    if (is.null(GISTIC) || !is.data.frame(GISTIC) || nrow(GISTIC)==0 ||
        !all(c("Type", "GISTIC_value", "From", "To") %in% colnames(GISTIC)))
        stop("Missing/invalid GISTIC levels definition data frame")

    # Gene names
    geneNames <- AnnotationDbi::select(
        org_obj,
        keytype="ENTREZID",
        keys=AnnotationDbi::keys(org_obj, keytype="ENTREZID"),
        columns=c("SYMBOL", "ENTREZID")
    )
    geneNames[,"ENTREZID"] <- as.character(geneNames[,"ENTREZID"])
    i <- (is.na(geneNames[,"SYMBOL"]) | geneNames[,"SYMBOL"]=="" |
              is.na(geneNames[,"ENTREZID"]) | geneNames[,"ENTREZID"]=="" |
              !grepl("^[0-9]+$", geneNames[,"ENTREZID"]) | geneNames[,"ENTREZID"]=="0")
    if (sum(i)>0) geneNames <- geneNames[!i,]

    # Extract genes from the transcript database & make sure that there is a symbol for all of them
    geneLocii <- GenomicFeatures::genes(tx_obj)
    i         <- which(as.character(GenomicRanges::mcols(geneLocii)$gene_id) %in% as.character(geneNames$ENTREZID))
    if (length(i) == 0) stop("No valid gene")
    geneLocii <- geneLocii[i,,drop=FALSE]

    # Find chromosome lengths
    chrlen <- BSgenome::seqinfo(bs_obj)
    GenomeInfoDb::seqlevelsStyle(chrlen) <- "NCBI"
    chrlen <- GenomeInfoDb::seqlengths(chrlen)
    assertthat::assert_that(all(ratios$Chromosome %in% names(chrlen)))

    # Find and read Control-FREEC output file
    bed <- collapse_ratio(ratios, collapse_col="CopyNumber", chrlen=chrlen)
    bed <- bed[bed[,"CopyNumber"]!=(-1),]
    if (nrow(bed) == 0)
        stop("No valid segment")
    bed$chrom <- paste("chr", bed$chrom, sep="")
    bed$GENES <- "."

    # Allele number to GISTIC
    tmp <- rep(as.integer(NA), nrow(bed))
    for (i in 1:nrow(GISTIC)) {
        j <- (GISTIC$From[i] <= bed$CopyNumber) & (bed$CopyNumber <= GISTIC$To[i])
        if (any(j)) tmp[j] <- GISTIC$GISTIC_value[i]
    }
    bed$CopyNumber <- tmp

    bedGR <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE, ignore.strand=TRUE)

    # Find overlaps with genes
    ovl <- GenomicRanges::findOverlaps(bedGR, geneLocii, ignore.strand=TRUE)
    ovSplit <- S4Vectors::split(S4Vectors::subjectHits(ovl), S4Vectors::queryHits(ovl))
    ovSplit <- sapply(ovSplit, function(j) {
        x <- as.character(geneLocii$gene_id[j])
        x <- x[!is.na(x)]
        paste(sort(unique(x)), collapse=";")
    })
    bed$GENES[as.integer(names(ovSplit))] <- ovSplit

    # subset to CNVs with annotated genes:
    genes <- bed %>%
        dplyr::mutate(GENES=strsplit(GENES, ";")) %>%
        tidyr::unnest(GENES)
    genes <- genes %>%
        dplyr::filter(GENES!=".") %>%
        dplyr::group_by(GENES) %>%
        dplyr::filter(CopyNumber==min(CopyNumber)) %>%
        dplyr::filter(dplyr::row_number()==1) %>%
        dplyr::ungroup()
    genes <- genes %>%
        dplyr::left_join(geneNames %>% dplyr::mutate(ENTREZID=as.character(ENTREZID)),
                         by=c("GENES"="ENTREZID"))
    genes <- genes %>%
        dplyr::select(Hugo_Symbol=SYMBOL, Entrez_Gene_Id=GENES, CopyNumber)

    # Store the calls
    Calls <- genes %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, CopyNumber)

    # Find and read Control-FREEC output file
    bed <- collapse_ratio(ratios, collapse_col="MedianRatio", chrlen=chrlen)
    bed <- bed[bed[,"MedianRatio"]!=(-1),]
    if (nrow(bed) == 0)
        stop("No valid segment")
    bed$chrom <- paste("chr", bed$chrom, sep="")
    bed$GENES <- "."
    bed$MedianRatio[bed$MedianRatio==0] <- small_value

    bedGR <- GenomicRanges::makeGRangesFromDataFrame(bed, keep.extra.columns=TRUE, ignore.strand=TRUE)

    # Find overlaps with genes
    ovl <- GenomicRanges::findOverlaps(bedGR, geneLocii, ignore.strand=TRUE)
    ovSplit <- S4Vectors::split(S4Vectors::subjectHits(ovl), S4Vectors::queryHits(ovl))
    ovSplit <- sapply(ovSplit, function(j) {
        x <- as.character(geneLocii$gene_id[j])
        x <- x[!is.na(x)]
        paste(sort(unique(x)), collapse=";")
    })
    bed$GENES[as.integer(names(ovSplit))] <- ovSplit

    # subset to CNVs with annotated genes:
    genes <- bed %>%
        dplyr::mutate(GENES=strsplit(GENES, ";")) %>%
        tidyr::unnest(GENES)
    genes <- genes %>%
        dplyr::filter(GENES!=".") %>%
        dplyr::group_by(GENES) %>%
        dplyr::summarise(MedianRatio=median(log2(MedianRatio))) %>%
        dplyr::ungroup()
    genes <- genes %>%
        dplyr::left_join(geneNames %>% dplyr::mutate(ENTREZID=as.character(ENTREZID)),
                         by=c("GENES"="ENTREZID"))
    genes <- genes %>%
        dplyr::select(Hugo_Symbol=SYMBOL, Entrez_Gene_Id=GENES, MedianRatio)

    # Store the log_2 medians
    Log_2 <- genes %>% dplyr::select(Hugo_Symbol, Entrez_Gene_Id, MedianRatio)

    # Segment file
    Segments <- bed %>%
        dplyr::mutate(chrom=sub("^chr", "", chrom),
                      num.mark=bins,
                      MedianRatio=log2(MedianRatio)) %>%
        dplyr::select(chrom, loc.start=start, loc.end=end, num.mark, seg.mean=MedianRatio) %>%
        dplyr::filter(!is.na(seg.mean) & !is.infinite(seg.mean))

    # CNVkit files
    CNVkit_cnv <- bed %>%
        dplyr::mutate(chrom=sub("^chr", "", chrom)) %>%
        dplyr::rename(chromosome=chrom, log2=MedianRatio) %>%
        dplyr::mutate(start=start-1, end=end, gene=".") %>%
        dplyr::select(chromosome, start, end, gene, log2) %>%
        dplyr::mutate(log2=log2(log2))

    CNVkit_cnr <- do.call("rbind",
                          lapply(split(ratios, ratios$Chromosome),
                                 function(x) {
                                     x <- x[order(x$Start),]
                                     x$end <- c(x$Start[-1]-1, chrlen[x$Chromosome[1]])
                                     x
                                 }))
    CNVkit_cnr <- CNVkit_cnr %>%
        dplyr::mutate(Start=Start-1, gene=".") %>%
        dplyr::select(chromosome=Chromosome, start=Start, end, gene, log2=Ratio) %>%
        dplyr::filter(log2!=(-1)) %>%
        dplyr::mutate(log2=replace(log2, log2==0, small_value)) %>%
        dplyr::mutate(log2=log2(log2)) %>%
        dplyr::arrange(chromosome, start) %>%
        as.data.frame()

    list(Calls=Calls, Log_2=Log_2, Segments=Segments, CNVkit_cnv=CNVkit_cnv, CNVkit_cnr=CNVkit_cnr)
}

#' Create files for pipeline & cnvkit plotting function from Control-FREEC ratio.txt file
#'
#' @param ratios_fn Control-FREEC ratio.txt file name
#' @param sample_name Sample name
#' @param log2_fn gene-base log2 CNV value (output file name)
#' @param call_fn GISTIC-like gene-based call (output file name)
#' @param segments_fn log2 CNV values for segment (output file name)
#' @param cns_fn cnvkit cns output file name
#' @param cnr_fn cnvkit cnr output file name
#' @param org_obj Bioconductor object defining ids mappings in an organism (default human)
#' @param tx_obj Bioconductor TxDb object defining gene locations (default human)
#' @param bs_obj Bioconductor BSgenome object defining the genome sequence(default human)
#'
#' @return TRUE invisibly. The side effects are writing the five files
#' @export
#'
#' @examples
control_freec_write_files <- function(sample_name,
                                      ratios_fn=paste("freec.", sample_name, "_ratio.txt", sep=""),
                                      log2_fn=paste("freec.", sample_name, "_gene_log2.txt", sep=""),
                                      call_fn=paste("freec.", sample_name, "_gene_call.txt", sep=""),
                                      segments_fn=paste("freec.", sample_name, "_segments.txt", sep=""),
                                      cns_fn=paste("freec.", sample_name, ".cns", sep=""),
                                      cnr_fn=paste("freec.", sample_name, ".cnr", sep=""),
                                      org_obj=org.Hs.eg.db::org.Hs.eg.db,
                                      tx_obj=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      bs_obj=BSgenome.Hsapiens.1000genomes.hs37d5::hs37d5) {
    if (missing(sample_name) || !is.character(sample_name) || length(sample_name)!=1 ||
        is.na(sample_name) || sample_name=="")
        stop("Missing or invalide sample name")
    if (is.null(ratios_fn) || !is.character(ratios_fn) || length(ratios_fn)!=1 ||
        is.na(ratios_fn) || ratios_fn=="" || file.access(ratios_fn, mode=4)!=0)
        stop("Missing, invalid or unreadable Control-FREEC ratios.txt file")
    if (is.null(log2_fn) || !is.character(log2_fn) || length(log2_fn)!=1 ||
        is.na(log2_fn) || log2_fn=="")
        stop("Missing or invalid gene-base log2 values output file name")
    if (is.null(call_fn) || !is.character(call_fn) || length(call_fn)!=1 ||
        is.na(call_fn) || call_fn=="")
        stop("Missing or invalid gene-based calls output file name")
    if (is.null(segments_fn) || !is.character(segments_fn) || length(segments_fn)!=1 ||
        is.na(segments_fn) || segments_fn=="")
        stop("Missing or invalid segments output file name")
    if (is.null(cns_fn) || !is.character(cns_fn) || length(cns_fn)!=1 ||
        is.na(cns_fn) || cns_fn=="")
        stop("Missing or invalid .cns output file name")
    if (is.null(cnr_fn) || !is.character(cnr_fn) || length(cnr_fn)!=1 ||
        is.na(cnr_fn) || cnr_fn=="")
        stop("Missing or invalid .cnr output file name")

    # Read Control-FREEC ratio.bed file
    ratios <- read.table(ratios_fn, sep="\t", header=1, stringsAsFactors=FALSE)

    # Compute everything
    tmp <- freec_to_calls_log2_segments(ratios, org_obj=org_obj, tx_obj=tx_obj, bs_obj=bs_obj)

    # Rename columns & add sample column to segments
    colnames(tmp$Log_2)[3] <- sample_name
    colnames(tmp$Calls)[3] <- sample_name
    tmp$Segments <- data.frame(ID=sample_name, tmp$Segments, check.names=FALSE, stringsAsFactors=FALSE)

    # Write output files
    write.table(tmp$Log_2, file=log2_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    write.table(tmp$Calls, file=call_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    write.table(tmp$Segments, file=segments_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    write.table(tmp$CNVkit_cnv, file=cns_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)
    write.table(tmp$CNVkit_cnr, file=cnr_fn, sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

    # Return nothing
    invisible(TRUE)
}

#' Read & merge sample files into cBioPortal gene-based data table.
#'
#' Used for expression tables, log2 CNA, pseudo-Gistic table.
#'
#' @param fns named list or vector of sample files, with sample ids as names
#' @param mappings mapping table between feature ids (must contain ENSEMBL, SYMBOL & ENTREZ_ID columns)
#' @param type type of data to merge
#' @param args named list of additional arguments required to process the data before exporting
#' @param remove_feature_version remove the trailing ".<number>" in the feature id (useful for GENCODE)
#'
#' @example
merge_tables <- function(fns, mappings, type=c("log2", "gistic", "segment", "expression"), args=list(), remove_feature_version=TRUE) {
    type <- match.arg(type)

    if (missing(fns) || is.null(fns))
        stop("Missing list of sample files")
    if (is.list(fns)) {
        stopifnot(all(lengths(fns) == 1))
        fns <- unlist(fns)
    }
    if (!is.character(fns) || any(is.na(fns)) || any(fns==""))
        stop("Missing at least one sample filename")
    if (is.null(names(fns)) || any(is.na(names(fns))) || any(names(fns)==""))
        stop("Missing at least one sample id")

    if (type == "segment") return(merge_segments(fns))

    if (missing(mappings) || is.null(mappings) || !is.character(mappings) || length(mappings)!=1 || is.na(mappings) || mappings=="")
        stop("Missing or illegal gene id mappings filename")
    mappings <- get_id_mappings(mappings, verbose=TRUE)

    if (type == "expression") {
        stopifnot(all(c("pipeline_id", "tx_obj") %in% names(args)))
        counts <- read_sample_files(fns, featureCol="V1", valueCol="V2", header=FALSE)
        tmp <- compute_rpkm(counts, tx_obj=args$tx_obj, verbose=TRUE)
        method <- "sum"
    }

    if (type == "gistic") {
        stopifnot(all(c("pipeline_id", "amplification") %in% names(args)))
        # Copy numbers (in "cn" column) are transformed into (pseudo-) gistic codes:
        # 0: Deep deletion, 1: heterozygous deletion, 2: copy number neutral, 3: gain, 4: amplification
        # In https://doi.org/10.1038/s41586-022-04738-6, the amplification is defined as 
        # copy number greater or equal to 9 copies (args$amplification = 9).
        args$amplification <- as.numeric(args$amplification)
        cn <- round(read_sample_files(fns, args$pipeline_id, "cn"))
        cn[args$amplification<=cn] <- args$amplification
        cn[2<cn & cn<args$amplification] <- 3
        cn[cn==args$amplification] <- 4
        # The gistic codes are changed so that copy number neutral alterations are at 0.
        # This is allowing finer selection of representative genes when there is no 1-1 mapping between pipeline id & gene symbols.
        # This way, one could select as representative the gene with more extreme copy number change, using "maxabs".
        # At the moment, the option describe above is not the selected one.
        # The representative is currently selected based on the highest copy number over all samples (sum of cn is highest).
        tmp <- cn-2
        method <- "max"
    }

    if (type == "log2") {
        stopifnot(all(c("pipeline_id") %in% names(args)))
        tmp <- read_sample_files(fns, args$pipeline_id, "log2")
        method <- "max"
    }

    if (remove_feature_version) rownames(tmp) <- sub("\\.[0-9]+$", "", rownames(tmp))

    rslt <- map_feature_id(tmp, mappings, from=args$pipeline_id, to="SYMBOL", method=method)

    rslt <- data.frame(
        Hugo_Symbol=rownames(rslt),
        Entrez_Gene_Id=mappings$ENTREZ_ID[match(rownames(rslt), mappings$SYMBOL)],
        rslt,
        row.names=NULL, check.names=FALSE
    )
    # Fix potentially missing ENTREZ_ID
    rslt[is.na(rslt[["Entrez_Gene_id"]]),"Entrez_Gene_Id"] <- 0

    rslt
}

#' Aggregate segments
#'
#' NOTE: this should not be done in R, python should be fine, but
#'       the csv package cannot handle long tokens produced by cnvkit in the gene column.
#' NOTE: the column renaming should be done by the copy number step.
#'
#' @param fns named list or named vector of sample filenames. The names are the sample IDs
#'
#' @example
merge_segments <- function(fns) {
    if (missing(fns) || is.null(fns))
        stop("Missing or illegal segmentation filenames")
    if (is.list(fns)) {
        stopifnot(all(lengths(fns) == 1))
        fns <- unlist(fns)
    }
    if (!is.character(fns) || any(is.na(fns)) || any(fns==""))
        stop("Missing file name")
    if (is.null(names(fns)))
        stop("Missing sample names")

    col_names <- c("chrom", "loc.start", "loc.end", "num.marks", "seg.mean")
    tbl <- NULL
    for (sample_id in names(fns)) {
        cat("Loading", fns[sample_id], "for sample", sample_id, "\n")
        tmp <- read.table(fns[sample_id], sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
        stopifnot(all(col_names %in% colnames(tmp)))
        tmp <- tmp[,col_names]
        tmp$ID <- sample_id
        tbl <- rbind(tbl, tmp)
    }

    return(tbl)
}

#' Aggregates counts and computes RPKM
#'
#' The computation of RPKM is done using DESeq2::fpkm in robust mode (after library size estimation).
#' The gene lengths are estimated using the total span of all its exons, which will lead
#' to gene length over-estimation in most cases.
#'
#' @param counts counts matrix or data.frame, feature ids in row names, sample ids in col names
#' @param tx_obj Bioconductor feature description object or its name (human by default, either TxDb or EnsDb), or GTF/GFF features file
#' @param verbose show progress of the computation
#'
#' @return the RPKM data frame
#' @export
#'
#' @examples
compute_rpkm <- function(counts, tx_obj=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene, verbose=FALSE) {
    if (missing(counts) || is.null(counts) || !(is.data.frame(counts) || is.matrix(counts)) || nrow(counts)<1 || ncol(counts)<1)
        stop("Missing or illegal counts data")
    tx_obj <- get_tx_object(tx_obj)

    counts <- counts[rowSums(!is.na(counts))>0,,drop=FALSE]
    if (nrow(counts) == 0)
        stop("No feature with counts")
    counts <- counts[,colSums(is.na(counts))==0,drop=FALSE]
    if (ncol(counts) == 0)
        stop("No sample without missing value")

    # Append the gene locii to the counts in dds object
    if (verbose) cat("Create DESeq2 object ... ")
    genes <- GenomicFeatures::exonsBy(tx_obj, "gene")
    genes <- genes[names(genes) %in% rownames(counts)]
    counts <- counts[names(genes),,drop=FALSE]
    donors <- data.frame(Donor=colnames(counts), stringsAsFactors=FALSE)
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData=donors, design=as.formula("~ 1"), rowRanges=genes)
    if (verbose) cat("Done\n")

    # Compute FPKMs
    if (verbose) cat("Compute FPKM ... ")
    dds <- DESeq2::estimateSizeFactors(dds)
    FPKM <- DESeq2::fpkm(dds, robust=TRUE)
    if (verbose) cat("Done\n")

    FPKM
}

require(magrittr)

#' Aggregates counts and computes RPKM
#'
#' The computation of RPKM is done using DESeq2::fpkm in robust mode (after library size estimation).
#' The gene lengths are estimated using the total span of all its exons, which will lead
#' to gene length over-estimation in most cases.
#'
#' @param df data frame describing expression counts & CNV gene calls files
#' @param count_input_format gene expression counts file format (only "feature_counts" supported)
#' @param tx_obj Bioconductor feature description object (human by default, either TxDb or EnsDb)
#' @param org_obj Bioconductor organism object (human by default)
#' @param pipeline_id type of ids used by the pipeline ("ENSEMBL" by default)
#' @param symbol_id Hugo_Symbol id name in org_obj ("SYMBOL" by default)
#' @param verbose show progress of the computation
#'
#' @return the RPKM data frame
#' @export
#'
#' @examples
compute_rpkm <- function(df,
                         count_input_format="feature_counts",
                         tx_obj=TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
                         org_obj=org.Hs.eg.db::org.Hs.eg.db, pipeline_id="ENSEMBL", symbol_id="SYMBOL",
                         verbose=FALSE) {
    if (missing(df) || is.null(df) ||
        !is.data.frame(df) || nrow(df)<1 ||
        !all(c("ID", "count_filename", "gene_call_filename") %in% colnames(df)))
        stop("Missing or illegal directory path")
    if (any(is.na(df$ID) | df$ID=="") || any(duplicated(df$ID)))
        stop("Illegal ID column")
    if (is.null(count_input_format) ||
        !is.character(count_input_format) || length(count_input_format)!=1 || is.na(count_input_format))
        stop("Missing or illegal file format for gene expression counts")
    count_input_format <- match.arg(count_input_format)
    if (is.null(tx_obj) || !(is(tx_obj, "TxDb") || is(tx_obj, "EnsDb")))
        stop("Missing or illegal features description object")
    if (is.null(org_obj) || !is(org_obj, "OrgDb"))
        stop("Missing or illegal organism description object")
    if (is.null(pipeline_id) || !is.character(pipeline_id) || length(pipeline_id)!=1 || is.na(pipeline_id) ||
        pipeline_id=="" || !(pipeline_id %in% AnnotationDbi::columns(org_obj)))
        stop("Missing, illegal or unknown pipeline id type")
    if (is.null(symbol_id) || !is.character(symbol_id) || length(symbol_id)!=1 || is.na(symbol_id) ||
        symbol_id=="" || !(symbol_id %in% AnnotationDbi::columns(org_obj)))
        stop("Missing, illegal or unknown Hugo_Symbol id type")

    # Identify samples with expression files
    has_expr <- !is.na(df$count_filename) & df$count_filename!="" & file.access(df$count_filename, mode=4)==0
    names(has_expr) <- as.character(df$ID)
    if (!any(has_expr))
        stop("No valid gene expression count file")

    # Read counts counts
    if (verbose) cat("Read expression counts ... ")
    counts <- switch(
        count_input_format,
        feature_counts=read_featurecounts_files(df$count_filename[has_expr])
    )
    colnames(counts) <- df$ID[has_expr]
    if (verbose) cat("Done\n")

    # Remove samples with no expression counts (all NA)
    has_expr[colnames(counts)] <- has_expr[colnames(counts)] & (colSums(!is.na(counts))>0)
    if (!any(has_expr))
        stop("No gene expression sample with valid expression counts")
    counts <- counts[,colSums(!is.na(counts))>0,drop=FALSE]

    # Remap counts to ENTREZ gene ids
    if (verbose) cat("Convert gene ids from ENSEMBL to ENTREZ ... ")
    rownames(counts) <- sub("\\.[0-9]+$", "", rownames(counts))
    counts <- pipeline_to_entrez(counts, org_obj=org_obj, pipeline_id=pipeline_id)
    if (verbose) cat("Done\n")

    # Append the gene locii to the counts in dds object
    if (verbose) cat("Create DESeq2 object ... ")
    genes <- GenomicFeatures::exonsBy(tx_obj, "gene")
    genes <- genes[names(genes) %in% rownames(counts)]
    counts <- counts[names(genes),]
    donors <- data.frame(Donor=colnames(counts), stringsAsFactors=FALSE)
    dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData=donors, design=as.formula("~ 1"), rowRanges=genes)
    if (verbose) cat("Done\n")

    # Compute FPKMs
    if (verbose) cat("Compute FPKM ... ")
    dds <- DESeq2::estimateSizeFactors(dds)
    FPKM <- DESeq2::fpkm(dds, robust=TRUE)
    if (verbose) cat("Done\n")

    # Make data frame with Hugo_Symbol, Entrez_Gene_Id and z-scores columns
    if (verbose) cat("Add gene symbols ... ")
    gene_names <- AnnotationDbi::select(org_obj,
                                        keys=rownames(FPKM),
                                        keytype="ENTREZID",
                                        columns=c(symbol_id, "ENTREZID"))
    if (is.null(gene_names) || !is.data.frame(gene_names) || nrow(gene_names)<1)
        stop("Entrez_Gene_Id ids not found")
    gene_names <- gene_names[!is.na(gene_names[[symbol_id]]) & !is.na(gene_names[["ENTREZID"]]),]
    if (nrow(gene_names)==0)
        stop("No valid mappings")
    gene_names <- gene_names[,c(symbol_id, "ENTREZID")]
    colnames(gene_names) <- c("Hugo_Symbol", "Entrez_Gene_Id")

    FPKM <- FPKM[rownames(FPKM) %in% as.character(gene_names$Entrez_Gene_Id),,drop=FALSE]
    i <- match(rownames(FPKM), as.character(gene_names$Entrez_Gene_Id))
    FPKM <- data.frame(gene_names[i,], FPKM, stringsAsFactors=FALSE, check.names=FALSE)
    rownames(FPKM) <- NULL
    if (verbose) cat("Done\n")

    FPKM
}

#' Aggregates counts and creates z-scores file
#'
#' Pipeline-implementation agnostic, location of files defined in the input data.frame
#'
#' @param df data frame describing expression counts & CNV gene calls files
#' @param count_input_format gene expression counts file format (only "feature_counts" supported)
#' @param gene_call_input_format CNV gene calls file format (only "cubi_pipeline" supported)
#' @param org_obj Bioconductor organism object (human by default)
#' @param pipeline_id type of ids used by the pipeline ("ENSEMBL" by default)
#' @param symbol_id Hugo_Symbol id name in org_obj ("SYMBOL" by default)
#' @param min_nb_diploid minimum number of diploid samples in order to compute z-scores (default 10)
#' @param verbose show progress of the computation
#' @param fast use vst instead of rlog to compute normalised expression
#'
#' @return the z_score data frame
#' @export
#'
#' @examples
compute_z_scores <- function(df,
                             count_input_format="feature_counts",
                             gene_call_input_format="cubi_pipeline",
                             org_obj=org.Hs.eg.db::org.Hs.eg.db, pipeline_id="ENSEMBL", symbol_id="SYMBOL",
                             min_nb_diploid=10,
                             verbose=FALSE, fast=FALSE) {
    if (missing(df) || is.null(df) ||
        !is.data.frame(df) || nrow(df)<1 ||
        !all(c("ID", "count_filename", "gene_call_filename") %in% colnames(df)))
        stop("Missing or illegal directory path")
    if (any(is.na(df$ID) | df$ID=="") || any(duplicated(df$ID)))
        stop("Illegal ID column")
    if (is.null(count_input_format) ||
        !is.character(count_input_format) || length(count_input_format)!=1 || is.na(count_input_format))
        stop("Missing or illegal file format for gene expression counts")
    count_input_format <- match.arg(count_input_format)
    if (is.null(gene_call_input_format) ||
        !is.character(gene_call_input_format) || length(gene_call_input_format)!=1 || is.na(gene_call_input_format))
        stop("Missing or illegal file format for CNV gene calls")
    gene_call_input_format <- match.arg(gene_call_input_format)
    if (is.null(org_obj) || !is(org_obj, "OrgDb"))
        stop("Missing or illegal organism description object")
    if (is.null(pipeline_id) || !is.character(pipeline_id) || length(pipeline_id)!=1 || is.na(pipeline_id) ||
        pipeline_id=="" || !(pipeline_id %in% AnnotationDbi::columns(org_obj)))
        stop("Missing, illegal or unknown pipeline id type")
    if (is.null(symbol_id) || !is.character(symbol_id) || length(symbol_id)!=1 || is.na(symbol_id) ||
        symbol_id=="" || !(symbol_id %in% AnnotationDbi::columns(org_obj)))
        stop("Missing, illegal or unknown Hugo_Symbol id type")

    # Identify samples with expression files
    has_expr <- !is.na(df$count_filename) & df$count_filename!="" & file.access(df$count_filename, mode=4)==0
    names(has_expr) <- as.character(df$ID)
    if (!any(has_expr))
        stop("No valid gene expression count file")

    # Read counts counts
    if (verbose) cat("Read expression counts ... ")
    counts <- switch(
        count_input_format,
        feature_counts=read_featurecounts_files(df$count_filename[has_expr])
    )
    colnames(counts) <- df$ID[has_expr]
    if (verbose) cat("Done\n")

    # Remove samples with no expression counts (all NA)
    has_expr[colnames(counts)] <- has_expr[colnames(counts)] & (colSums(!is.na(counts))>0)
    if (!any(has_expr))
        stop("No gene expression sample with valid expression counts")
    counts <- counts[,colSums(!is.na(counts))>0,drop=FALSE]

    # Remap counts to ENTREZ gene ids
    if (verbose) cat("Convert gene ids from ENSEMBL to ENTREZ ... ")
    rownames(counts) <- sub("\\.[0-9]+$", "", rownames(counts))
    counts <- pipeline_to_entrez(counts, org_obj=org_obj, pipeline_id=pipeline_id)
    if (verbose) cat("Done\n")

    # Normalise counts
    if (verbose) cat("Normalise expression ... ")
    has_expr[colnames(counts)] <- has_expr[colnames(counts)] & (colSums(is.na(counts))==0)
    counts <- counts[,colSums(is.na(counts))==0,drop=FALSE] # Required for DESeq2::DESeq
    if (!any(has_expr))
        stop("No valid gene expression sample for normalisation")
    expr <- normalise_counts(counts, verbose=verbose, fast=fast)
    if (verbose) cat("Done\n")

    # Identify files with CNV gene calls and
    has_cnv <- !is.na(df$gene_call_filename) & df$gene_call_filename!="" & file.access(df$gene_call_filename, mode=0)==0
    has_both <- has_expr & has_cnv

    cnv_calls <- NULL
    if (any(has_both)) {
        if (verbose) cat("Read CNV gene calls ... ")
        cnv_calls <- switch(
            gene_call_input_format,
            cubi_pipeline=process_cnv_genes(df[has_both,] %>% dplyr::rename(filename=gene_call_filename))
        )
        rownames(cnv_calls) <- as.character(cnv_calls$Entrez_Gene_Id)
        cnv_calls <- cnv_calls %>% dplyr::select(-Hugo_Symbol, -Entrez_Gene_Id) %>% as.matrix()
        colnames(cnv_calls) <- df$ID[has_both]
        if (verbose) cat("Done\n")
    }

    # Compute z-scores (when no CNV calls are present, use all samples)
    if (verbose) cat("Compute expression z-scores ... ")
    z <- z_scores(expr, cnv_calls, min_nb_diploid=min_nb_diploid)
    if (verbose) cat("Done\n")

    # Gene symbols from ENTREZ ids
    if (verbose) cat("Add gene symbols ... ")
    if (any(has_both)) {
        gene_names <- process_cnv_genes(df[has_both,] %>% dplyr::rename(filename=gene_call_filename)) %>%
            dplyr::select(Hugo_Symbol, Entrez_Gene_Id)
    } else {
        if (!all(c(symbol_id, "ENTREZID") %in% AnnotationDbi::columns(org_obj)))
            stop("Missing/illegal organism database (OrgDb object)")

        gene_names <- AnnotationDbi::select(org_obj,
                                            keys=rownames(z),
                                            keytype="ENTREZID",
                                            columns=c(symbol_id, "ENTREZID"))
        if (is.null(gene_names) || !is.data.frame(gene_names) || nrow(gene_names)<1)
            stop("Entrez_Gene_Id ids not found")
        gene_names <- gene_names[!is.na(gene_names[[symbol_id]]) & !is.na(gene_names[["ENTREZID"]]),]
        if (nrow(gene_names)==0)
            stop("No valid mappings")
        gene_names <- gene_names[,c(symbol_id, "ENTREZID")]
        colnames(gene_names) <- c("Hugo_Symbol", "Entrez_Gene_Id")
    }
    if (verbose) cat("Done\n")

    # Make data frame with Hugo_Symbol, Entrez_Gene_Id and z-scores columns
    z <- z[rownames(z) %in% as.character(gene_names$Entrez_Gene_Id),,drop=FALSE]
    i <- match(rownames(z), as.character(gene_names$Entrez_Gene_Id))
    z <- data.frame(gene_names[i,], z, stringsAsFactors=FALSE, check.names=FALSE)
    rownames(z) <- NULL

    z
}

#' Compute gene expression z-scores from diploid samples
#'
#' @param norm_expr normalised expression matrix
#' @param cnv_calls matrix of calls (0 for diploid)
#' @param min_nb_diploid minimum number of diploid samples in order to compute z-scores (default 10)
#'
#' @return a matrix of z scores
#' @export
#'
#' @examples
z_scores <- function(norm_expr, cnv_call=NULL, min_nb_diploid=10) {
    if (missing(norm_expr) || is.null(norm_expr) ||
        !is.matrix(norm_expr) || nrow(norm_expr)<1 || ncol(norm_expr)<1 ||
        !is.numeric(norm_expr) || any(is.na(norm_expr)) ||
        is.null(rownames(norm_expr)) || any(is.na(rownames(norm_expr))) || any(rownames(norm_expr)=="") ||
        is.null(colnames(norm_expr)) || any(is.na(colnames(norm_expr))) || any(colnames(norm_expr)==""))
        stop("Missing/illegal normalised expression matrix")
    if (is.null(min_nb_diploid) || !is.numeric(min_nb_diploid) || length(min_nb_diploid)!=1 ||
        is.na(min_nb_diploid) || min_nb_diploid<0)
        stop("Missing/illegal minimum number of diploid samples")
    if (min_nb_diploid<3) min_nb_diploid <- 3

    filter <- matrix(TRUE, nrow=nrow(norm_expr), ncol=ncol(norm_expr), dimnames=dimnames(norm_expr))

    if (!is.null(cnv_call)) {
        if (!is.matrix(cnv_call) || nrow(cnv_call)<1 || ncol(cnv_call)<1 ||
            is.null(rownames(cnv_call)) || any(is.na(rownames(cnv_call))) || any(rownames(cnv_call)=="") ||
            is.null(colnames(cnv_call)) || any(is.na(colnames(cnv_call))) || any(colnames(cnv_call)==""))
            stop("Illegal gene-based copy-number matrix")
        cnv_call <- !is.na(cnv_call) & cnv_call==0
        cnv_call <- cnv_call[rownames(cnv_call) %in% rownames(norm_expr),colnames(cnv_call) %in% colnames(norm_expr),drop=FALSE]
        if (nrow(cnv_call)>0 & ncol(cnv_call)>0)
            filter[rownames(cnv_call),colnames(cnv_call)] <- cnv_call
    }

    tmp <- norm_expr
    tmp[!filter] <- as.numeric(NA)
    N <- rowSums(!is.na(tmp))
    m <- rowMeans(tmp, na.rm=TRUE)
    s <- apply(tmp, 1, sd, na.rm=TRUE)

    s[s==0 | N<min_nb_diploid] <- as.numeric(NA)
    (norm_expr - m)/s
}

#' Read feature_count output files
#'
#' @param fns feature counts output file names
#'
#' @return a numeric matrix of counts, with gene ids as row names
#' @export
#'
#' @examples
read_featurecounts_files <- function(fns) {
    if (missing(fns) || !is.character(fns) || length(fns)==0 ||
        sum(!is.na(fns) & fns!="")==0)
        stop("Missing/illegal featurecounts file names")

    tmp <- lapply(fns, function(x) {
        if (is.na(x) || x=="" || !file.exists(x) || file.access(x, mode=4)!=0)
            return(data.frame(Geneid=as.character(), Counts=as.integer(), stringsAsFactors=FALSE))
        read.table(x, sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
    })
    stopifnot(all(sapply(tmp, function(x) all(x[,1]==tmp[[1]][,1]))))
    stopifnot(all(sapply(tmp, function(x) all(x[,2]==tmp[[1]][,2]))))
    stopifnot(all(sapply(tmp, function(x) all(x[,3]==tmp[[1]][,3]))))
    stopifnot(all(sapply(tmp, function(x) all(x[,4]==tmp[[1]][,4]))))
    stopifnot(all(sapply(tmp, function(x) all(x[,5]==tmp[[1]][,5]))))
    stopifnot(all(sapply(tmp, function(x) all(x[,6]==tmp[[1]][,6]))))

    ids <- unique(unlist(lapply(tmp, function(x) x$Geneid)))
    ids <- ids[!is.na(ids) & as.character(ids)!=""]
    if( length(ids) == 0)
        stop("No readable non-empty featurecounts file")

    counts <- matrix(as.numeric(NA), nrow=length(ids), ncol=length(tmp),
                     dimnames=list(ids, fns))
    for (i in 1:length(tmp)) {
        if (!is.null(tmp[[i]]) && is.data.frame(tmp[[i]]) && nrow(tmp[[i]]) > 0) {
            tmp[[i]] <- tmp[[i]][tmp[[i]]$Geneid %in% ids,,drop=FALSE]
            if (nrow(tmp[[i]])<1 || any(is.na(tmp[[i]]$Geneid)) || any(duplicated(tmp[[i]]$Geneid)))
                next
            counts[tmp[[i]]$Geneid,i] <- tmp[[i]][,ncol(tmp[[i]])]
        }
    }

    counts
}

#' Remap ids for gene expression count matrix
#'
#' Given a gene expression count matrix with ids stored as row names, the function maps these ids to ENTREZ ids,
#' and outputs a gene expression count matrix with ENTREZ ids as rownames.
#'
#' When several original ids map to the same ENTREZ id, the output value will be the median of each input entry.
#'
#' @param pip_counts Gene expression count matrix, with gene ids as row names
#' @param org_obj Bioconductor organism description object containing mappings between original & ENTREZ gene ids (default: human)
#' @param pipeline_id type of ids stored in the pip_counts row names ("ENSEMBL" by default)
#' @param FUNC which value to use in presence of multiple pipeline genes for one ENTREZ gene (default: max)
#'
#' @return gene expression count matrix for ENTREZ genes
#' @export
#'
#' @examples
pipeline_to_entrez <- function(pip_counts, org_obj=org.Hs.eg.db::org.Hs.eg.db, pipeline_id="ENSEMBL", FUNC=max) {
    if (missing(pip_counts) || is.null(pip_counts) ||
        !is.matrix(pip_counts) || nrow(pip_counts)<1 || ncol(pip_counts)<1 ||
        !is.numeric(pip_counts) ||
        is.null(rownames(pip_counts)) || any(is.na(rownames(pip_counts))) || any(rownames(pip_counts)==""))
        stop("Missing/illegal counts matrix")
    if (is.null(pipeline_id) || !is.character(pipeline_id) || length(pipeline_id)!=1 || is.na(pipeline_id) ||
        pipeline_id=="")
        stop("Missing, illegal or unknown pipeline id type")
    if (pipeline_id == "ENTREZID") return(pip_counts)
    if (is.null(org_obj) || !is(org_obj, "OrgDb") ||
        !all(c(pipeline_id, "ENTREZID") %in% AnnotationDbi::columns(org_obj)))
        stop("Missing/illegal organism database (OrgDb object)")

    ids <- AnnotationDbi::select(org_obj,
                                 keys=rownames(pip_counts),
                                 keytype=pipeline_id,
                                 columns=c(pipeline_id, "ENTREZID"))
    if (is.null(ids) || !is.data.frame(ids) || nrow(ids)<1)
        stop("Pipeline ids not found")
    ids <- ids[!is.na(ids[,pipeline_id]) & !is.na(ids[,"ENTREZID"]),]
    if (nrow(ids)==0)
        stop("No valid mappings")
    ids[,pipeline_id] <- as.character(ids[,pipeline_id])
    ids[,"ENTREZID"] <- as.character(ids[,"ENTREZID"])

    i <- split(ids[,pipeline_id], ids[,"ENTREZID"])
    i1 <- i[lengths(i)==1]
    i2 <- i[lengths(i)>1]

    ent_counts <- c()
    if (length(i1)>0) {
        i1 <- unlist(i1)
        tmp <- pip_counts[i1,]
        rownames(tmp) <- names(i1)
        ent_counts <- rbind(ent_counts, tmp)
    }
    if (length(i2)>0) {
        tmp <- t(sapply(i2, function(j) apply(pip_counts[j,,drop=FALSE], 2, FUNC)))
        ent_counts <- rbind(ent_counts, tmp)
    }

    round(ent_counts)
}

#' Normalise gene expression counts on regularised log scale
#'
#' @param counts gene expression count matrix (must be integer)
#' @param verbose DESeq2 normalisation progress messages
#' @param fast use vst instead of rlog to compute normalised expression
#'
#' @return matrix of normalised expression on regularised log scale
#' @export
#'
#' @examples
normalise_counts <- function(counts, verbose=FALSE, fast=FALSE) {
    if (missing(counts) || is.null(counts) ||
        !is.matrix(counts) || nrow(counts)<1 || ncol(counts)<1 ||
        !is.numeric(counts) || any(is.na(counts)) || any(counts<0) ||
        is.null(rownames(counts)) || any(is.na(rownames(counts))) || any(rownames(counts)=="") ||
        is.null(colnames(counts)) || any(is.na(colnames(counts))) || any(colnames(counts)==""))
        stop("Missing/illegal counts matrix")

    donors <- data.frame(Donor=colnames(counts), stringsAsFactors=FALSE)

    normalised <- DESeq2::DESeqDataSetFromMatrix(
        countData=counts, colData=donors, design=as.formula("~ 1"), tidy=FALSE
    )
    normalised <- DESeq2::DESeq(normalised, betaPrior=TRUE, quiet=!verbose)
    if (fast) {
        normalised <- DESeq2::vst(normalised, blind=TRUE)
    } else {
        normalised <- DESeq2::rlog(normalised, blind=TRUE)
    }
    SummarizedExperiment::assay(normalised)
}

#' Merge gene-based CNV results
#'
#' @param df gene calls file names
#'
#' @return data.frame with one row per gene and one column per donor.
#'     The first two columns are Hugo_Symbol (HGNC gene name) & Entrez_Gene_Id.
#'     Only genes with positive ENTREZ ids are returned.
#' @export
#'
#' @examples
process_cnv_genes <- function(df) {
    if (missing(df) || is.null(df) || !is.data.frame(df) || nrow(df)==0)
        stop("Missing or illegal CNV files data frame")

    df <- df %>%
        dplyr::filter(!is.na(filename) & file.access(filename, mode=4)==0)

    unprocessed <- c()
    rslt <- NULL
    for (i in 1:nrow(df)) {
        x <- read.table(df$filename[i], sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
        x <- x[rowSums(is.na(x))==0,]

        # Workaround a bug in the first versions of the copywriter target_seq_cnv_calling step
        if (colnames(x)[1] == "Hugo_Gene_Symbol") colnames(x)[1] <- "Hugo_Symbol"

        if (ncol(x) != 3 || !all(colnames(x)[1:2] == c("Hugo_Symbol", "Entrez_Gene_Id"))) {
            unprocessed <- c(unprocessed, df$ID[i])
            next
        }

        if (is.null(rslt)) {
            rslt <- x[,1:2]
            ref_hash <- paste(rslt$Hugo_Symbol, rslt$Entrez_Gene_Id, sep="\t")
        }

        x_hash <- paste(x$Hugo_Symbol, x$Entrez_Gene_Id, sep="\t")

        rslt <- rslt[ref_hash %in% x_hash,]
        ref_hash <- ref_hash[ref_hash %in% x_hash]
        x <- x[match(ref_hash, x_hash),]

        rslt <- cbind(rslt, x[,3])
        colnames(rslt)[ncol(rslt)] <- df$ID[i]
    }

    if (length(unprocessed) > 0)
        warning("Gene CNA data could not be processed for donors ", paste(unprocessed, collapse=", "))

    if (all(grepl("^ *[+-]?[0-9]+ *$", rslt$Entrez_Gene_Id))) {
        rslt$Entrez_Gene_Id <- as.integer(rslt$Entrez_Gene_Id)
        if ("Entrez_Gene_Id" %in% colnames(rslt))
            rslt <- rslt %>% dplyr::filter(Entrez_Gene_Id>0)
    } else {
        warning("Non-numerical Entrez_Gene_Id records")
    }

    if (any(duplicated(rslt$Entrez_Gene_Id))) {
        warning("Duplicated ENTREZ genes ",
                paste(unique(rslt$Entrez_Gene_Id[duplicated(rslt$Entrez_Gene_Id)]), collapse=", "))
        rslt <- rslt[!duplicated(rslt$Entrez_Gene_Id),]
    }

    rslt
}


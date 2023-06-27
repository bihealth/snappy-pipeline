require(magrittr)

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
merge_tables <- function(fns, mappings, type=c("log2", "gistic", "expression"), args=list(), remove_feature_version=TRUE) {
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
        cn <- read_sample_files(fns, args$pipeline_id, "cn")
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
    if (type == "gistic") rslt <- rslt+2

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
    if ("chromosome" %in% colnames(segments)) colnames(segments)[colnames(segments)=="chromosome"] <- "chrom"
    stopifnot(all(c("cn", "log2", "chrom", "start", "end") %in% colnames(segments)))

    seg <- GenomicRanges::makeGRangesFromDataFrame(segments, seqnames.field="chrom", start.field="start", end.field="end", strand="*", keep.extra.columns=FALSE)

    i <- GenomicRanges::findOverlaps(genes, seg)
    i <- split(S4Vectors::subjectHits(i), S4Vectors::queryHits(i))
    i <- unlist(i[lengths(i) == 1])

    tbl <- data.frame(FeatureID=as.character(genes$gene_id), cn=NA, log2=NA)
    colnames(tbl)[1] <- pipeline_id
    tbl[as.numeric(names(i)),c("cn", "log2")] <- segments[i,c("cn", "log2")]

    tbl
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

    col_names <- c("chromosome"="chrom", "start"="loc.start", "end"="loc.end", "x"="num.mark", "log2"="seg.mean")
    tbl <- NULL
    for (sample_id in names(fns)) {
        tmp <- read.table(fns[i], sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE)
        stopifnot(all(names(col_names) %in% colnames(tmp)))
        tmp <- tmp[,names(col_names)]
        colnames(tmp) <- col_names    
        tmp <- cbind(ID=rep(sample_id, nrow(tmp)), tmp)
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
    counts <- counts[names(genes),]
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

#' Extracts the feature names & one column of values from tab-delimited sample files
#' and return a data.frame with the values for each sample in the columns, and 
#' the feature names as row name.
#'
#' @param fns named list or named vector of sample filenames. The names are the sample IDs
#' @param featureCol the name of the column containing the feature names in the sample files
#' @param valueCol the name of the column containing the values to be extracted from each sample file.
#'
#' @example
read_sample_files <- function(fns, featureCol, valueCol, header=TRUE) {
    if (missing(fns) || is.null(fns))
        stop("Missing list of files")
    if (is.list(fns)) {
        stopifnot(all(lengths(fns) == 1))
        fns <- unlist(fns)
    }
    if (!is.character(fns) || any(is.na(fns)) || any(fns==""))
        stop("Missing file name")
    if (is.null(names(fns)))
        stop("Missing sample names")
    if (missing(featureCol) || is.null(featureCol) || !is.character(featureCol) || length(featureCol)!=1 || is.na(featureCol) || featureCol=="")
        stop("Missing or illegal feature column name")
    if (missing(valueCol) || is.null(valueCol) || !is.character(valueCol) || length(valueCol)!=1 || is.na(valueCol) || valueCol=="")
        stop("Missing or illegal column name for data extraction")

    tbl <- NULL
    for (i in names(fns)) {
        tmp <- read.table(fns[i], sep="\t", header=header, stringsAsFactors=FALSE, check.names=FALSE)
        stopifnot(all(c(featureCol, valueCol) %in% colnames(tmp)))
        tmp <- tmp[,c(featureCol, valueCol),drop=FALSE]
        colnames(tmp)[2] <- i
        if (is.null(tbl)) {
            tbl <- tmp
        } else {
            tbl <- tbl %>% dplyr::full_join(tmp, by=featureCol)
        }
    }
    tbl <- tbl[!grepl("_PAR_Y$", tbl[[featureCol]]),,drop=FALSE]
    if (any(duplicated(tbl[[featureCol]])))
        stop("Duplicated feature ids")
    rownames(tbl) <- tbl[[featureCol]]
    tbl[,colnames(tbl)!=featureCol,drop=FALSE]
}


#' Remap data to different feature ID (typically from ENSEMBL to SYMBOL or ENTREZ_ID).
#'
#' When the mappings between ids is not 1-1, the following rules are applied:
#' many-many relationships are removed
#' when pipeline id maps to multiple symbols, the data is copied to all symbols
#' when multiple pipeline ids map to a single symbol, the the method is taken into account:
#' when method is sum, the data for the target id is the sum of all "source" id for every sample,
#' when method is max, the data of the "source" id which is highest over all samples is selected as representative for the "target" id.
#'
#' @param mat the matrix or data.frame containing the data (expression, CNA, ...)
#' @param mappings the data frame of feature id mappings
#' @param from the "source" type of feature ids, must be in the mapping table column names
#' @param to the "target" type of feature ids, must be in the mapping table column names
#' @param method aggregation method when many "source" ids are mapped to the same "target" id
#'
#' @example
map_feature_id <- function(mat, mappings, from, to, method=c("sum", "max", "maxabs")) {
    if (missing(mat) || is.null(mat) || !(is.data.frame(mat) || is.matrix(mat)) || nrow(mat)<1 || ncol(mat)<1)
        stop("Missing or illegal data")
    method <- match.arg(method)
    if (missing(from) || is.null(from) || !is.character(from) || length(from)!=1 || is.na(from) || from=="")
        stop("Missing or illegal feature id type (origin)")
    if (missing(to) || is.null(to) || !is.character(to) || length(to)!=1 || is.na(to) || to=="")
        stop("Missing or illegal feature id type (target)")
    if (from==to) return(mat)
    if (missing(mappings) || is.null(mappings) || !is.data.frame(mappings) || nrow(mappings)<1 || ncol(mappings)<2)
        stop("Missing or illegal feature id mapping table")
    if (!all(c(from, to) %in% colnames(mappings)))
        stop("Feature id type ", from, " or ", to, " not in mappings table")

    mappings <- mappings[,c(from, to)] %>% dplyr::distinct()
    mappings <- mappings[mappings[[from]] %in% rownames(mat),,drop=FALSE]
    
    n <- sum(rownames(mat) %in% mappings[[from]])
    if (n == 0)
        stop("No common feature id between the mapping table and the data")
    if (n/nrow(mat) < 0.9)
        warning("The proportion of feature id that can be mapped is low at ", n/nrow(mat))
    if (n < nrow(mat))
        warning(nrow(mat) - n, " features cannot be mapped, removed from data")
    mat <- mat[rownames(mat) %in% mappings[[from]],,drop=FALSE]
    mappings <- mappings[mappings[[from]] %in% rownames(mat),,drop=FALSE]
    
    # Resolve 1-1, 1-many & many-1 mappings between ENSEMBL & HGNC symbols
    mappings <- map_between_two_ids(mappings[[from]], mappings[[to]], haveIndex=FALSE, verbose=TRUE, keepBlanks=FALSE)
    rslt <- NULL
    if (any(mappings[["type"]] == "1-1")) {
        i <- which(mappings[["type"]] == "1-1")
        tmp <- mat[mappings[["x"]][i],,drop=FALSE]
        rownames(tmp) <- mappings[["y"]][i]
        rslt <- rbind(rslt, tmp)
    }
    if (any(mappings[["type"]] == "1-many")) {
        i <- which(mappings[["type"]] == "1-many")
        tmp <- mat[mappings[["x"]][i],,drop=FALSE]
        rownames(tmp) <- mappings[["y"]][i]
        rslt <- rbind(rslt, tmp)
    }
    if (any(mappings[["type"]] == "many-1")) {
        i <- which(mappings[["type"]] == "many-1")
        i <- split(mappings[["x"]][i], mappings[["y"]][i])
        tmp <- switch(method,
            sum=sapply(i, function(j) colSums(mat[j,,drop=FALSE], na.rm=TRUE)),
            max=sapply(i, function(j) { x <- mat[j,,drop=FALSE] ; x <- x[order(-rowSums(x, na.rm=TRUE)),,drop=FALSE] ; x[1,] }),
            maxabs=sapply(i, function(j) { x <- mat[j,,drop=FALSE] ; x <- x[order(-rowSums(abs(x), na.rm=TRUE)),,drop=FALSE] ; x[1,] })
        )
        if (!is.matrix(x)) tmp <- matrix(tmp, ncol=1, dimnames=list(names(i), colnames(mat)))
        rslt <- rbind(rslt, tmp)
    }

    rslt
}

#' Gets a TxDb/EnsDb object
#'
#' The TxDb/EnsDb object can be obtained from different inputs
#' - An existing TxDb or EnsDb object, in this case the function is only a pass-through,
#' - The name of TxDb/EnsDb package, possibly qualified with the object name in the environment, or
#' - The name of a GTF/GFF file, or
#' - An existing GRanges object
#'
#' @param tx_obj TxDb/EnsDb, GRanges, filename or package/object name
#'
#' @return TxDb/EnsDb object
#' @export
#'
#' @examples
get_tx_object <- function(tx_obj, verbose=FALSE) {
    if (missing(tx_obj) || is.null(tx_obj))
        stop("Missing features description object")
    if (is.character(tx_obj) && length(tx_obj)==1 && !is.na(tx_obj) && tx_obj!="") {
        if (file.exists(tx_obj)) {
            if (verbose) cat("Taking features from GTF/GTF file ", tx_obj)
            tx_obj <- GenomicFeatures::makeTxDbFromGFF(tx_obj)
        } else {
            if (grepl("::", tx_obj)) {
                tx_obj <- getFromNamespace(sub("^.+::", "", tx_obj), ns=sub("::.+$", "", tx_obj))
            } else {
                tx_obj <- getFromNamespace(tx_obj, ns=tx_obj)
            }
        }
    } else {
        if (is(tx_obj, "GRanges")) tx_obj <- GenomicFeatures::makeTxDbFromGRanges(tx_obj)
    }
    if (!(is(tx_obj, "TxDb") || is(tx_obj, "EnsDb")))
        stop("Illegal features description object")
    tx_obj
}

#' Gets a gene id mapping table
#'
#' The mapping table contains gene identifiers from HGNC (SYMBOL), ENSEMBL & NCBI (ENTREZ_ID).
#' The mapping table has no missing identifier, and no duplicate rows.
#' All ENTREZ_ID are numbers coded as characters.
#'
#' Several type of inputs are possible:
#' - An existing orgDb object, 
#' - The name of an orgDb package, possibly qualified with the object name in the environment, or
#' - A tab-delimited file, containined the columns "SYMBOL", "ENSEMBL" and "ENTREZ_ID". 
#'   To enable the possibility to use genome-nexus files, a file with column names
#'   "hgnc_symbol", "ensembl_canonical_gene", "entrez_gene_id" is also allowed.
#'
#' @param org_obj orgDb, mapping id filename or package/object name
#'
#' @return gene id mapping table with columns SYMBOL, ENSEMBL & ENTREZ_ID
#' @export
#'
#' @examples
get_id_mappings <- function(org_obj, verbose=FALSE) {
    if (missing(org_obj) || is.null(org_obj))
        stop("Missing organism description object")
    if (is(org_obj, "OrgDb")) {
        if (verbose) cat("Gene id mappings taken from Bioconductor object")
        id_mappings <- AnnotationDbi::select(
            org_obj, keytype=pipeline_id, columns=c("ENSEMBL", "SYMBOL", "ENTREZ_ID"),
            keys=AnnotationDbi::keys(org_obj, keytype=pipeline_id)
        )
    } else {
        if (!is.character(org_obj) || length(org_obj)!=1 || is.na(org_obj) || org_obj=="")
            stop("Illegal organism description")
        if (file.exists(org_obj)) {
            if (verbose) cat("Gene id mappings taken from file ", org_obj)
            id_mappings <- read.table(org_obj, sep="\t", header=1, stringsAsFactors=FALSE, check.names=FALSE, quote="", comment="")
            if (all(c("hgnc_symbol", "ensembl_canonical_gene", "entrez_gene_id") %in% colnames(id_mappings)))
                id_mappings <- id_mappings %>%
                    dplyr::select(ENSEMBL=ensembl_canonical_gene, SYMBOL=hgnc_symbol, ENTREZ_ID=entrez_gene_id)
            stopifnot(all(c("ENSEMBL", "SYMBOL", "ENTREZ_ID") %in% colnames(id_mappings)))
            id_mappings <- id_mappings[,c("ENSEMBL", "SYMBOL", "ENTREZ_ID")]
        } else {
            if (grepl("::", org_obj)) {
                org_obj <- getFromNamespace(sub("^.+::", "", org_obj), ns=sub("::.+$", "", org_obj))
            } else {
                org_obj <- getFromNamespace(org_obj, ns=org_obj)
            }
            id_mappings <- AnnotationDbi::select(
                org_obj, keytype=pipeline_id, columns=c("ENSEMBL", "SYMBOL", "ENTREZ_ID"),
                keys=AnnotationDbi::keys(org_obj, keytype=pipeline_id)
            )
        }
    }
    id_mappings <- id_mappings %>%
        dplyr::select(ENSEMBL, SYMBOL, ENTREZ_ID) %>%
        dplyr::mutate(ENTREZ_ID=as.character(ENTREZ_ID)) %>%
        dplyr::filter(!is.na(ENSEMBL) & ENSEMBL!="") %>%
        dplyr::filter(!is.na(SYMBOL) & SYMBOL!="") %>%
        dplyr::filter(!is.na(ENTREZ_ID) & ENTREZ_ID!="" & grepl("^[0-9]+$", ENTREZ_ID)) %>%
        dplyr::distinct()

    id_mappings
}

#' Create table of mappings between two sets of ids
#'
#' Creates a complete mapping between two sets of ids (x & y) of the same length.
#' The mappings are given a type (1-1, 1-many, many-1 & many-many), and
#' a group identifier, so that the mappings of all y ids corresponding to
#' a single x id in a 1-many relationship are given the same group number.
#' Missing values and duplicate entries in the input are ignored.
#'
#' The function creates a data frame containing all mappings between ids in
#' set x and ids in set y. These mappings are given a type (1-1, 1-many,
#' many-1 and many-many), and the original mappings are completed to include
#' all mappings between an x and a y which are connected by some path.
#' For example, if the x and y inputs are \code{x = c("x1", "x2", "x2")} and
#' \code{y = c("y1", "y1", "y2")}, which describe that id x1 is mapped to y1,
#' x2 is mapped to y1 and y2, then in the output, x1 will be mapped to y2.
#' A single group id will be assigned to the 4 connections (x1 to y1, x1 to y2,
#' x2 to y1 and x2 to y2).
#' If \code{haveIndex = TRUE}, the function will also create a column which
#' contains the indices of the original mappings in the input vectors. In the
#' example above, the index column will be 1 for x1 to y1, 2 for x2 to y1,
#' 3 for x2 to y2 and NA for x1 to y2, as the latter is not present in the
#' input.
#' Note that the index column is a column of lists, as it will report all
#' occurences of the original mapping, which might appear several times.
#' The missing mappings (which have NAs in either x or y input vectors) are
#' discarded.
#'
#' @param x vector of ids
#' @param y vector of ids
#' @param haveIndex return index of mapping in original vectors
#' @param verbose output progress messages
#' @param keepBlanks allow the empty string as valid ID
#'
#' @return data.frame with columns x (the first id), y (the second id),
#'     type (the mapping type, 1-1, 1-many, many-1 or many-many),
#'     group (an integer grouping id pairs), and optionally index, a
#'     pointer to the mapping's position in the initial input. If the
#'     mapping is not in the original input (i.e. it has been created
#'     because of other connections between x and y ids), then the value
#'     of the index is NA.
#' @export
#'
#' @examples{
#' map_between_two_ids(c("x1", "x2", "x2", "x3", "x4", "x5", "x6", "x6"),
#'                     c("y1", "y1", "y2", "y3", "y4", "y4", "y5", "y6"))
#' }
map_between_two_ids <- function( x, y, haveIndex=TRUE, verbose=TRUE, keepBlanks=FALSE ) {
    # Cleanup input
    if( missing(x) || missing(y) ||
        is.null( x ) || is.null( y ) || length( x ) < 1 ||
        length( x ) != length( y ) )
        stop( "Missing arguments" )
    if( verbose ) cat( "Starting grouping", length( x ), "elements\n" )
    df <- unique( cbind( x=as.character( x ), y=as.character( y ) ) )
    rownames( df ) <- NULL
    df <- df[rowSums( is.na( df ) )==0,]
    if( ! keepBlanks ) df <- df[rowSums( df == "" )==0,]
    if( verbose ) cat( "Number of complete pairs:", nrow( df ), "\n" )
    if( nrow( df ) == 0 ) return( NULL )

    # Split by x and y columns
    ys_from_x <- split( df[,"y"], df[,"x"] )
    xs_from_y <- split( df[,"x"], df[,"y"] )
    one_y_from_x <- unlist( ys_from_x[lengths(ys_from_x) == 1] )
    many_ys_from_x <- ys_from_x[lengths(ys_from_x) > 1]
    one_x_from_y <- unlist( xs_from_y[sapply( xs_from_y, length )==1] )
    many_xs_from_y <- xs_from_y[sapply( xs_from_y, length )>1]

    igroup <- 0
    # 1-1 mappings: those in one_y_from_x that also are in one_x_from_y
    rslt <- data.frame( group=as.numeric( c() ),
                        x=as.character( c() ), y=as.character( c() ),
                        type=as.character(),
                        stringsAsFactors=FALSE )
    if( length( one_y_from_x ) > 0 && length( one_x_from_y ) > 0 ) {
        tmp <- one_y_from_x[names( one_y_from_x ) %in% one_x_from_y]
        if (length(tmp)>0) {
            tmp <- data.frame( group=igroup+(1:length( tmp )),
                               x=names( tmp ), y=tmp,
                               type="1-1",
                               stringsAsFactors=FALSE )
            rownames( tmp ) <- NULL
            igroup <- nrow(tmp)
            rslt <- rbind(rslt, tmp)

            # Remove 1-1 mappings from single splits
            one_y_from_x <- one_y_from_x[!(one_y_from_x %in% tmp$y)]
            one_x_from_y <- one_x_from_y[!(one_x_from_y %in% tmp$x)]
        }
        if( verbose ) cat( "Number of 1-1 pairs:", nrow( tmp ), "\n" )
    }

    # many-1 mappings: each y has many xs (only xs which map to a single y must be considered)
    if( length( one_y_from_x ) > 0 && length( many_xs_from_y ) > 0 ) {
        many_to_1_x <- names( one_y_from_x )
        many_to_one <- many_xs_from_y[sapply( many_xs_from_y, function( y ) all( y %in% many_to_1_x ) )]
        n <- 0
        if (length(many_to_one) > 0) {
            tmp <- data.frame( group=igroup + rep( (1:length( many_to_one )), lengths( many_to_one ) ),
                               x=unlist(many_to_one), y=rep(names(many_to_one), lengths(many_to_one)),
                               type="many-1",
                               stringsAsFactors=FALSE )
            rownames( tmp ) <- NULL
            rslt <- rbind( rslt, tmp )
            n <- nrow(tmp)
            igroup <- igroup + length(many_to_one)

            # Remove many-1 mappings from single splits
            many_xs_from_y <- many_xs_from_y[!(names( many_xs_from_y ) %in% names( many_to_one ))]
            one_y_from_x <- one_y_from_x[!(names(one_y_from_x) %in% many_to_1_x)]
        }
        if( verbose ) cat( "Number of many-1 pairs:", n, "\n" )
    }

    # 1-many mappings: each x has many ys (only ys that map to a single xmust be considered)
    if( length( one_x_from_y ) > 0 && length( many_ys_from_x ) > 0 ) {
        one_to_many_y <- names( one_x_from_y )
        one_to_many <- many_ys_from_x[sapply( many_ys_from_x, function( x ) all( x %in% one_to_many_y ) )]
        n <- 0
        if (length(one_to_many) > 0) {
            tmp <- data.frame( group=igroup + rep( (1:length( one_to_many )), lengths( one_to_many ) ),
                               x=rep(names(one_to_many), lengths(one_to_many)), y=unlist(one_to_many),
                               type="1-many",
                               stringsAsFactors=FALSE )
            rownames( tmp ) <- NULL
            rslt <- rbind( rslt, tmp )
            n <- nrow(tmp)
            igroup <- igroup + length(one_to_many)

            many_ys_from_x <- many_ys_from_x[!(names( many_ys_from_x ) %in% names( one_to_many ))]
            one_x_from_y <- one_x_from_y[!(names(one_x_from_y) %in% one_to_many_y)]
        }
        if( verbose ) cat( "Number of 1-many pairs:", n, "\n" )
    }

    # many-many mappings
    if( length( many_ys_from_x )>0  && length(many_xs_from_y)>0) {
        # Create a matrix with all remaining connections
        m <- rbind( data.frame(x=rep(names(many_ys_from_x), lengths(many_ys_from_x)), y=unlist(many_ys_from_x),
                               stringsAsFactors=FALSE),
                    data.frame(x=unlist(many_xs_from_y), y=rep(names(many_xs_from_y), lengths(many_xs_from_y)),
                               stringsAsFactors=FALSE),
                    data.frame(x=one_x_from_y, y=names(one_x_from_y), stringsAsFactors=FALSE),
                    data.frame(x=names(one_y_from_x), y=one_y_from_x, stringsAsFactors=FALSE)
        )
        m <- as.matrix( unique( m ) )
        if( verbose ) cat( "Number of many-many pairs:", nrow( m ), "\n" )
        if( verbose ) cat( "Finding groups...\n" )

        # Get the connected components for this graph
        m[,1] <- paste( "x", m[,1], sep="\t" )
        m[,2] <- paste( "y", m[,2], sep="\t" )
        g     <- igraph::graph_from_edgelist( m, directed=FALSE )
        cl    <- igraph::components( g, mode="strong" )
        i     <- split( names( cl$membership ), cl$membership )
        assertthat::assert_that(min(lengths(i))>3)
        if( verbose ) cat( "Number of many-many groups:", length( i ), "\n" )

        # Transform so that each group in i has 2 lists which contain the x & y items
        i     <- lapply( i, function( ii ) matrix( unlist( stringr::str_split( ii, "\t" ) ), nrow=2 ) )
        i     <- lapply( i, function( ii ) split( ii[2,], ii[1,] ) )

        # Build the dataframe for all, filling it with all combinations of x & y
        n    <- sum( sapply( i, function( ii ) length( ii$x )*length( ii$y ) ) )
        tmp  <- data.frame( group=rep( NA, n ), x=rep( NA, n ), y=rep( NA, n ), type="many-many",
                            stringsAsFactors=FALSE)
        k    <- 0
        for( j in names( i ) ) {
            igroup <- igroup+1
            nx <- length( i[[j]]$x )
            ny <- length( i[[j]]$y )
            tmp[k+(1:(nx*ny)),"group"] <- igroup
            tmp[k+(1:(nx*ny)),"x"]     <- rep( i[[j]]$x, each=ny )
            tmp[k+(1:(nx*ny)),"y"]     <- rep( i[[j]]$y, nx )
            k <- k + nx*ny
        }
        rslt <- rbind( rslt, tmp )
    }

    # Show the mappings that were actually in the input
    if( haveIndex ) {
        i <- split( 1:length( x ), paste( as.character( x ), as.character( y ), sep="\t" ) )
        j <- match( paste( as.character( rslt$x ), as.character( rslt$y ), sep="\t" ), names( i ) )
        rslt$index <- c()
        rslt$index[!is.na( j )] <- i[j[!is.na( j )]]
        assertthat::assert_that(all(!is.na(rslt$index)[rslt$type!="many-many"]))
    }

    rslt
}

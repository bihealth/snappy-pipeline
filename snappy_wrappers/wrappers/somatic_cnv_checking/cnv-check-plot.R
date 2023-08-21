cn_to_call <- function(cn, amplification=7) {
    Call <- factor(rep("Neutral", length(cn)), levels=c("Deletion", "Loss", "Neutral", "Gain", "Amplification"))
    Call[!is.na(cn) & cn==0] <- "Deletion"
    Call[!is.na(cn) & cn==1] <- "Loss"
    Call[!is.na(cn) & 2<cn & cn<amplification] <- "Gain"
    Call[!is.na(cn) & amplification<=cn] <- "Amplification"
    Call[is.na(cn)] <- NA
    Call
}

vcf_to_table <- function(vcf, sample, amplification=7) {
     vcf <- VariantAnnotation::readVcf(vcf)
     stopifnot(all(lengths(VariantAnnotation::alt(vcf)) == 2))
     
     x <- data.frame(
         CHROM=GenomicRanges::seqnames(SummarizedExperiment::rowRanges(vcf)),
         POS=GenomicRanges::start(SummarizedExperiment::rowRanges(vcf)),
         REF=VariantAnnotation::ref(vcf),
         ALT=unlist(lapply(VariantAnnotation::alt(vcf), function(a) a[1])),
         t(sapply(VariantAnnotation::geno(vcf)[["AD"]][,sample], function(a) a[1:2])),
         CN=VariantAnnotation::geno(vcf)[["CN"]][,sample],
         LFC=VariantAnnotation::geno(vcf)[["LFC"]][,sample]
     )
     colnames(x)[5:6] <- c("t_ref", "t_alt")
     x <- x |>
         dplyr::filter(!is.na(LFC) & !is.na(CN)) |>
         dplyr::mutate(BAF=t_alt/(t_ref+t_alt)) |>
         dplyr::mutate(Call=cn_to_call(CN))
    x
}

cnv_to_table <- function(cnv) {
    y <- read.table(cnv, sep="\t", header=0)
    colnames(y) <- c("CHROM", "start", "stop", "name", "LFC", "strand")
    y <- y |>
        dplyr::mutate(LFC=replace(LFC, strand == "-", -abs(LFC[strand=="-"]))) |>
        dplyr::mutate(start=start + 1)
    y
}

chromosome_lengths <- function(genome, use_fai=TRUE, use_dict=TRUE) {
    if (use_fai && file.exists(paste(genome, "fai", sep="."))) {
        genome <- read.table(paste(genome, "fai", sep="."), sep="\t", header=0)
        colnames(genome) <- c("CHROM", "Length", "index", "ncol", "nbyte")
    } else {
        if (use_dict && file.exists(paste(genome, "dict", sep="."))) {
            pattern <- "^SQ\tSN:([^ \t]+)\tLN:([0-9]+)\t.+$"
            genome <- grep(pattern, readLines(paste(genome, "dict", sep=".")), value=TRUE)
            genome <- data.frame(CHROM=sub(pattern, "\\1", genome), Length=as.numeric(sub(pattern, "\\2", genome)))
        } else {
            genome <- Biostrings::readDNAStringSet(genome)
            genome <- data.frame(CHROM=names(genome), Length=Biostrings::width(genome))
        }
    }
    genome <- genome |>
        dplyr::mutate(CHROM=sub(" .*", "", CHROM)) |>
        dplyr::mutate(Offset=cumsum(as.numeric(dplyr::lag(Length, default=0)))) |>
        dplyr::select(CHROM, Length, Offset)
    genome
}

default_call_colors <- c(
    "Deletion"="#E41A1C",
    "Loss"="#984EA3",
    "Neutral"="#4DAF4A",
    "Gain"="#80B1D3",
    "Amplification"="#377EB8"
)

plot_cnv <- function(x, scale=c("log2", "sqrt"), call_colors=default_call_colors) {
    scale <- match.arg(scale)
    if (scale == "sqrt") x_label <- "Segment Fold Change"
    else                 x_label <- "Segment Log2 Fold Change"
    if (scale == "sqrt") x <- x |> dplyr::mutate(LFC=2^LFC)
    p <- x |>
        dplyr::mutate(BAF=pmin(BAF, 1-BAF)) |>
        ggplot2::ggplot(ggplot2::aes(x=LFC, y=BAF, color=Call)) +
        ggplot2::geom_point() +
        ggplot2::labs(x=x_label, y="B-allele fraction") +
        ggplot2::scale_color_manual(values=call_colors) +
        ggplot2::theme_bw()
    if (scale == "sqrt") p <- p + ggplot2::scale_x_sqrt()
    p
}

plot_locus <- function(x, genome_lengths, call_colors=default_call_colors) {
    p <- genome_lengths |>
        ggplot2::ggplot(ggplot2::aes(xmin=Offset+1, xmax=Offset+Length, ymin=0, ymax=1, fill=as.character(n%%2))) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_manual(values=c("white", "grey")) +
        ggplot2::scale_x_continuous(breaks=genome_lengths$Offset+genome_lengths$Length/2, labels=genome_lengths$CHROM) +
        ggplot2::geom_point(data=x, mapping=ggplot2::aes(x=x, y=BAF, color=Call), inherit.aes=FALSE) +
        ggplot2::scale_color_manual(values=call_colors) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1)) +
        ggplot2::guides(fill="none") +
        ggplot2::labs(x="Genomic position", y="B-allele fraction")
    p
}

plot_segment <- function(y, genome_lengths, call_colors=default_call_colors) {
    p <- genome_lengths |>
        ggplot2::ggplot(ggplot2::aes(xmin=Offset+1, xmax=Offset+Length, ymin=0, ymax=0.5, fill=as.character(n%%2))) +
        ggplot2::geom_rect() +
        ggplot2::scale_fill_manual(values=c("white", "grey")) +
        ggplot2::scale_x_continuous(breaks=genome_lengths$Offset+genome_lengths$Length/2, labels=genome_lengths$CHROM)
    p <- p +
        ggplot2::geom_point(data=y, mapping=ggplot2::aes(x=(from+to)/2, y=Median, color=Call, size=Nvariants), inherit.aes=FALSE) +
        ggplot2::geom_segment(data=y, mapping=ggplot2::aes(x=(from+to)/2, y=`1st quartile`, xend=(from+to)/2, yend=`3rd quartile`, color=Call), inherit.aes=FALSE) +
        ggplot2::geom_segment(data=y, mapping=ggplot2::aes(x=from, y=Median, xend=to, yend=Median, color=Call), inherit.aes=FALSE) +
        ggplot2::scale_color_manual(values=call_colors)
    p <- p +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x=ggplot2::element_text(angle=90, vjust=0.5, hjust=1)) +
        ggplot2::guides(fill="none") +
        ggplot2::labs(x="Genomic position", y="B-allele fraction")
    p
}

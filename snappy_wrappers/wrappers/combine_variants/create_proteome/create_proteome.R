#' Clean GENCODE genome
#' 
#' @description Removes extra identifiers, contig lengths, md5 checksums, ...
#' and leaves only the main contig name (typically, chr1, ...)
#' Onyl contigs that are present in the feature annotation object
#' are retained
#' 
#' @param genome Biostrings::DNAStringSet object, typically read from fasta file downloaded from GENCODE (or other source)
#' @param tx TxDb object containing all feature annotations, typically created from a gff file downloaded from GENCODE
#' 
#' @return a Biostrings::DNAStringSet with only the subset of contigs described in the feature annotation, and matching contig names
#' 
clean_genome <- function(reference, tx) {
    names(reference) <- sub(" .*", "", names(reference))
    reference[names(reference) %in% names(GenomeInfoDb::seqinfo(tx))]
}

#' Mutates a genome reference
#' 
#' @description The reference is mutated as a list of loci (which MUST BE non-overlapping).
#' The mutated genome is returned, as well as the set of chromosome sequence
#' shifts which must be applied to the features genomic coordinates in order
#' to keep those genomic coordinates in sync with the mutated sequence.
#' 
#' @param reference Biostrings::DNAStringSet object containing the un-mutated, reference sequence
#' @param loci GenomicRanges::GRanges object containing the mutation loci, as well as the reference & alt allele sequences
#' 
#' @return a list with 2 components, the mutated genome sequence in "genome" (Biostrings::DNAStringSet), 
#' and a data frame with 3 columns: the contig name, the mutation position (end), and the shift.
#' 
#' @details the shift column represents the effect of the mutation on the genomic position
#' at the end of the mutation. It is the difference between the lengths of the reference and
#' alternative alleles. It is very important to use only non-overlapping variants, the routine will fail othrewise.
#' Note that only chromosomes that have been altered by mutations are returned.
#' 
.add_non_overlapping_mutations <- function(reference, loci) {
    genome <- reference
    altered_chr <- c()
    shifts <- NULL
    for (chr_name in names(reference)) {
        chr <- reference[chr_name]
        var_on_chr <- loci[GenomicRanges::seqnames(loci) == chr_name]
        if (length(var_on_chr) == 0) next
        altered_chr <- c(altered_chr, chr_name)
        ranges <- GenomicRanges::ranges(var_on_chr)
        chr <- Biostrings::replaceAt(chr, ranges, value=var_on_chr$ALT)
        shifts <- rbind(shifts, data.frame(
            seqnames=chr_name,
            pos=IRanges::end(ranges),
            delta=nchar(as.character(var_on_chr$ALT)) - IRanges::width(ranges)
        ))
        genome[chr_name] <- chr
    }
    list(genome=genome[altered_chr], shifts=shifts |> dplyr::filter(delta != 0))
}

#' Inserts (typically germline) variants in a reference genome
#' 
#' @description The variants are inserted in a reference genome.
#' The variants are first decomposed into non-overlapping sets, and variants in
#' each set are separately inserted in the genome. Therefore, the function returns
#' several distinct modified genomes.
#' The input variants must be normalised (no multiple alt alleles).
#' 
#' @param reference Biostrings::DNAStringSet object containing the un-mutated, reference sequence
#' @param variants VariantAnnotation::VCF object containing all variants (typically germline)
#' 
#' @return list of mutated genomes, one per set of non-overlapping variants.
#' Each mutated genome is itself a list of 2 components, the mutated genome sequence in 
#' "genome" (Biostrings::DNAStringSet), and a data frame with 3 columns: the contig name, 
#' the mutation position (end), and the shift.
#' 
#' @details The splitting of the variants into non-overlapping set is such that
#' for each overlapping variant, the first one (which starting position is closest to
#' the contig's 5'end) is put in the first set, and the other variant in the second set.
#' We have not yet encountered more than two overlapping variants.
#' 
#' @note These overlapping variants frequently result from heterozyous indels in close
#' proximity with other SNVs or indels. Those cases often confuse variant callers,
#' and the resulting calls are somtimes quite questionable. However, it is very
#' difficult (and dangerous) to attempt to merge or re-write the calls. Therefore,
#' the decision has been taken to consider the calls separately, and produce
#' distinct genome alterations.
#' 
insert_variants <- function(reference, variants) {
    stopifnot(all(lengths(VariantAnnotation::alt(variants)) == 1))
    loci <- SummarizedExperiment::rowRanges(variants)
    loci$ALT <- unlist(loci$ALT)

    genomes <- list()
    repeat {
        i <- as.data.frame(GenomicRanges::findOverlaps(loci, loci, ignore.strand=TRUE)) |>
            dplyr::filter(queryHits < subjectHits)
        postponed <- i$subjectHits
        if (length(postponed) > 0) {
            current_loci <- loci[-postponed]
            postponed_loci <- loci[postponed]
        } else {
            current_loci <- loci
            postponed_loci <- NULL
        }

        genomes <- c(genomes, list(.add_non_overlapping_mutations(reference, current_loci)))
        if (is.null(postponed_loci)) break
        loci <- postponed_loci
    }

    genomes
}

#' Shift feature genomic coordinates by differences due to indels.
#' 
#' @description For each contig separately, the features loci are updated by the cumulative
#' differences due to indels. Both the ends of the features are updated (but not necessarily
#' by the same value when an indel in within the feature).
#' 
#' @param loci the loci of all features as GenomicRanges::GRangeList object. The GenomicRanges::GRanges
#' members of the list usually describe the loci of CDS exons, grouped by transcripts.
#' @param shifts a data frame with individual shifts due to indels. It is usually produced
#' by "insert_variants".
#' 
#' @return the updated loci coordinates, grouped by transcript (GenomicRanges::GRangeList object)
#' 
.shift_coordinate <- function(loci, shifts) {
    x <- unlist(loci)
    x$tx_name <- rep(names(loci), lengths(loci))
    x <- GenomicRanges::split(x, GenomicRanges::seqnames(x))

    shifted_coords <- NULL
    for (chr_name in names(x)) {
        chr_shifts <- shifts |> 
            dplyr::filter(seqnames == chr_name) |>
            dplyr::arrange(pos) |>
            dplyr::mutate(diff=cumsum(delta))
        y = as.data.frame(x[[chr_name]])
        if (nrow(chr_shifts) > 0) {
            i <- findInterval(y$start, c(0, chr_shifts$pos))
            y$start = y$start + c(0, chr_shifts$diff)[i]
            i <- findInterval(y$end, c(0, chr_shifts$pos))
            y$end = y$end + c(0, chr_shifts$diff)[i]
        }
        shifted_coords <- rbind(shifted_coords, y)
    }
    shifted_coords <- GenomicRanges::makeGRangesFromDataFrame(shifted_coords, keep.extra.columns=TRUE)
    GenomicRanges::split(shifted_coords, shifted_coords$tx_name)
}

#' Create a proteome based on a single mutated genome
#' 
#' @description Based on the mutated genome and the indels shift, the function
#' creates a the complete list of mutated proteome. When the shifts are not
#' present, the function returns the amino-acid sequences based on original loci.
#' This can be used to create a "reference" proteome.
#' 
#' @param genome (generally mutated) contig sequences as Biostrings::DNAStringSet object
#' @param tx GenomicFeatures::TxDb database describing all genes, transcripts & exons.
#' This is generally obtained using GenomicFeatures::makeTxDbFromGFF on GENCODE gff file.
#' @param shifts when present, this must be a data frame describing the shifts due to 
#' mutation indels. When absent, the loci from the TxDb database are unchanged.
#' 
#' @return Set of proteins as Biostrings::AAStringSet
#' 
create_proteome <- function(genome, tx, shifts=NULL) {
    local_tx <- tx
    GenomeInfoDb::seqlevels(local_tx) <- names(genome)
    cds_loci <- GenomicFeatures::cdsBy(local_tx, by="tx", use.names=TRUE)
    if (!is.null(shifts)) cds_loci <- .shift_coordinate(cds_loci)
    cds_seq <- GenomicFeatures::extractTranscriptSeqs(genome, cds_loci)
    i <- Biostrings::vmatchPattern("N", cds_seq)
    cds_seq <- cds_seq[lengths(i) == 0]
    aa_seq <- Biostrings::translate(cds_seq)
    aa_seq[order(names(aa_seq))]
}

#' Create all proteins than can possibly generated from a set of mutated genomes & one annotation
#' 
#' @param genomes list of mutated genomes created by "insert_variants"
#' @param tx GenomicFeatures::TxDb database describing all genes, transcripts & exons.
#' @param reference reference (unmutated) genome sequence as Biostrings::DNAStringSet
#' 
#' @return Complete set of proteins as Biostrings::AAStringSet from all mutated genomes and reference if provided
#' 
donor_proteome <- function(genomes, tx, reference=NULL) {
    proteomes <- NULL
    for (genome in genomes) {
        proteome <- as.character(create_proteome(genome$genome, tx, shifts=genome$shifts))
        n <- sum(grepl("\\*.+", proteome))
        proteome <- sub("\\*.*", "", proteome)
        proteomes <- unique(c(proteomes, proteome))
    }
    if (!is.null(reference)) {
        proteome <- as.character(create_proteome(genome$genome, tx))
        proteome <- sub("\\*.*", "", proteome)
        proteomes <- unique(c(proteomes, proteome))
    }
    proteomes
}

#' Merge overlapping variants
#' 
#' @description Attempt to merge consecutive indels into a single sequence change.
#' Only deletion followed by inserts can be processed, and only when there is
#' one position in common. Other overlapping mutations are ignored, and won't be
#' used to generate a personal proteome.
#' 
#' @param variants the variants as GenomicRanges::GRanges object (generally created
#' using SummarizedExperiment::rowRanges() on a VariantAnnotation::VCF object)
#' 
#' @return a modified set of variants as a GenomicRanges::GRanges object
#' 
#' @note This function is not used, it has been tested on a limited number of cases,
#' but the variants produced appeared to be reasonable. However, the number of
#' discarded overlapping variants is too high, and the lack of validation of the
#' merged variants are reasons why this approach (fixing overlapping variants) is not used.
#' 
.fix_overlapping_variants <- function(variants) {
    iLoop <- 0
    repeat {
        iLoop <- iLoop + 1
        message("iLoop = ", iLoop)
        i <- GenomicRanges::findOverlaps(variants, variants, ignore.strand=TRUE)
        i <- as.data.frame(i[S4Vectors::queryHits(i) < S4Vectors::subjectHits(i)])
        if (nrow(i) == 0) break

        if (iLoop > 1) {
            l <- list(i=i)
            l$variants <- as.data.frame(variants)
            return(l)
        }

        added <- NULL
        removed <- NULL

        not_processed <- unique(c(
            i$queryHits[duplicated(i$queryHits)],
            i$subjectHits[duplicated(i$subjectHits)],
            intersect(i$queryHits, i$subjectHits)
        ))
        if (length(not_processed) > 0) {
            removed <- c(removed, not_processed)
            i <- i[-not_processed,]
        }

        if (nrow(i) > 0) {
            delins <- which(
                (Biostrings::width(variants$REF[i$queryHits]) > Biostrings::width(variants$ALT[i$queryHits])) &
                (Biostrings::width(variants$REF[i$subjectHits]) < Biostrings::width(variants$ALT[i$subjectHits]))
            )
            not_processed <- setdiff(1:nrow(i), delins)
            if (length(not_processed) > 0) {
                removed <- c(removed, not_processed)
                i <- i[-not_processed,]
            }

            if (nrow(i) > 0) {
                for (j in 1:nrow(i)) {
                    v1 <- variants[i$queryHits[j]]
                    v2 <- variants[i$subjectHits[j]]
                    alt <- substring(as.character(v2$ALT), GenomicRanges::end(v1) - GenomicRanges::start(v2) + 1)
                    alt <- paste(substring(as.character(v1$REF), 1, 1), substring(alt, 2), sep="")
                    v0 <- v1
                    v0$ALT <- Biostrings::DNAStringSet(alt)
                    v0$QUAL <- min(v1$QUAL, v2$QUAL)
                    names(v0) <- sprintf("%s:%s_%s/%s", GenomicRanges::seqnames(v0), GenomicRanges::start(v0), as.character(v0$REF), as.character(v0$ALT))
                    if (is.null(added)) added <- v0
                    else                added <- c(added, v0)
                }
                removed <- c(removed, i$queryHits, i$subjectHits)
            }
        }

        if (length(removed) > 0) variants <- variants[-unique(removed)]
        variants <- c(variants, added)
    }

    variants
}

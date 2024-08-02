# -*- coding: utf-8 -*-
"""Wrapper for running postprocessing of CNVetti-based segmentation."""

import os

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

shell.executable("/bin/bash")

rscript = os.path.join(
    os.path.dirname(os.path.realpath(__file__)), "../../../copywriter/call/snappy-run-copywriter.R"
)

awkscript = os.path.join(os.path.dirname(os.path.realpath(__file__)), "merge_segments.awk")

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# -------------------------------------------------------------------------------------------------
# Write segment and bin/target log2 values to tsv files to read into R
SAMPLE=$(bcftools query -l {snakemake.input.targets_bcf})

bcftools view -f .,PASS {snakemake.input.targets_bcf} \
   | bcftools query -f "$SAMPLE\t%ID\t%CHROM\t%POS\t%INFO/END[\t%SG2]\n" \
   > {snakemake.output.targets_segmented_txt}

bcftools view -f .,PASS {snakemake.input.targets_bcf} \
   | bcftools query -f "$SAMPLE\t%ID\t%CHROM\t%POS\t%INFO/END[\t%CV2]\n" \
   > {snakemake.output.targets_txt}

awk -v OFS='\t' -f {awkscript} {snakemake.output.targets_segmented_txt} \
   > {snakemake.output.segments_txt}

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
cat <<"EOF" > $TMPDIR/run_cnvetti_postprocessing.R
library(CGHcall)
source("{rscript}")

raws = read.delim("{snakemake.output.targets_txt}", header=FALSE, stringsAsFactors=FALSE)
segs = read.delim("{snakemake.output.targets_segmented_txt}", header=FALSE, stringsAsFactors=FALSE)

# careful, ID here references the sample ID, not the feature ID (Feature)
column_names = c("ID", "Feature", "chrom", "loc.start", "loc.end", "seg.mean")
colnames(raws) = column_names
colnames(segs) = column_names
rownames(raws) = raws$Feature
rownames(segs) = segs$Feature

sampleID = unlist(raws$ID[1])
donorID=paste(strsplit(sampleID, "-", fixed=TRUE)[[1]][1:2], collapse="-")

# the following code is ripped off from the work that Eric did for copywriter
# we might deduplicate this at a later stage
valCols <- which( !(colnames( raws ) %in% c( "chrom", "maploc", "loc.start", "loc.end", "iBin", "Feature", "ID", "sample" )) )
nas_to_drop = which( rowSums( is.na( raws[,valCols,drop=FALSE] ) ) == 0 )

df.features <- AnnotatedDataFrame( data.frame( Chromosome=as.numeric( raws$chrom[nas_to_drop] ),
                                               Start=raws$loc.start[nas_to_drop],
                                               End=raws$loc.end[nas_to_drop],
                                               row.names=raws$Feature[nas_to_drop],
                                               check.names=FALSE ),
                                   varMetadata=data.frame( labelDescription=c( "Chromosomal position",
                                                                               "Basepair position start",
                                                                               "Basepair position end" ),
                                                           row.names=c( "Chromosome", "Start", "End" ) ),
                                   dimLabels=c( "featureNames", "colnames" ) )

# dummy
df.pheno    <- AnnotatedDataFrame( data.frame( x=rep(0, 1),
                                               row.names=make.names("seg.mean") )[,-1],
                                   varMetadata=as.data.frame( matrix( NA, nrow=0, ncol=1,
                                                                      dimnames=list( NULL, "labelDescription" ) ) ) )

cghseg = new( "cghSeg", copynumber=raws[valCols], segmented=segs[valCols],
              featureData=df.features, phenoData=df.pheno, annotation="hg19" )
sgn <- postsegnormalize(cghseg)
cll <- CGHcall(sgn)
res <- ExpandCGHcall(cll, sgn)

pdf(paste("{snakemake.output.segments_txt}", ".pdf", sep=""))
plot(res, dotres=5)
dev.off()

# Overwrite one of the copywriter functions
geneCalls <- function(segData, res, binSize = 20000, type = c("call", "log2")) {{
    require(TxDb.Hsapiens.UCSC.hg19.knownGene)
    require(org.Hs.eg.db)
    type <- match.arg(type)
    geneNames <- AnnotationDbi::select(org.Hs.eg.db, columns = "SYMBOL",
        keys = keys(org.Hs.eg.db, keytype = "ENTREZID"), keytype = "ENTREZID")
    i <- which(is.na(geneNames$ENTREZID) | is.na(geneNames$SYMBOL) |
        geneNames$ENTREZID == "" | as.numeric(geneNames$ENTREZID) <=
        0 | geneNames$SYMBOL == "")
    if (length(i) > 0)
        geneNames <- geneNames[-i, , drop = FALSE]
    if (nrow(geneNames) == 0)
        stop("No genes")
    genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
    i <- which(mcols(genes)$gene_id %in% geneNames$ENTREZID)
    if (length(i) == 0)
        stop("No common genes")
    genes <- genes[i, , drop = FALSE]
    # change: set coltype from sample name to call or log2
    colnames(res) = type
    # change: match target df with calls
    seg <- to_GRanges(cbind(segData, call =unlist(calls(res))),
        binSize = binSize)
    sampleNames <- unique(segData$ID)
    rslt <- matrix(NA, nrow = length(genes), ncol = length(sampleNames),
        dimnames = list(mcols(genes)$gene_id, sampleNames))
    for (sample in sampleNames) {{
        s <- seg[mcols(seg)$ID == sample, ]
        i <- as.matrix(findOverlaps(genes, s))
        i <- split(i[, "subjectHits"], i[, "queryHits"])
        # take the first bin a gene matches to
        i = unlist(lapply(i, `[[`, 1))
        if (type == "call")
            rslt[as.numeric(names(i)), sample] <- mcols(s)$call[i]
        if (type == "log2")
            rslt[as.numeric(names(i)), sample] <- mcols(s)$seg.mean[i]
    }}
    i <- match(mcols(genes)$gene_id, geneNames$ENTREZID)
    geneNames <- geneNames[i, ]
    data.frame(Hugo_Symbol = geneNames$SYMBOL, Entrez_Gene_Id = geneNames$ENTREZID,
        rslt)
}}


# Output the calls
tmp <- geneCalls( segs, res, binSize=binSize, type="call" )
colnames( tmp )[3] <- donorID
write.table( tmp,
             file="{snakemake.output.gene_call_txt}",
             sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
tmp <- geneCalls( segs, res, binSize=binSize, type="log2" )
colnames( tmp )[3] <- donorID
write.table( tmp,
             file="{snakemake.output.gene_log2_txt}",
             sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

warnings()
sessionInfo()

save.image(paste("{snakemake.output.segments_txt}", ".RData", sep=""))
EOF

Rscript --vanilla $TMPDIR/run_cnvetti_postprocessing.R
md5sum {snakemake.output.segments_txt} > {snakemake.output.segments_txt}.md5
md5sum {snakemake.output.targets_txt} > {snakemake.output.targets_txt}.md5
md5sum {snakemake.output.targets_segmented_txt} > {snakemake.output.targets_segmented_txt}.md5
md5sum {snakemake.output.gene_call_txt} > {snakemake.output.gene_call_txt}.md5
md5sum {snakemake.output.gene_log2_txt} > {snakemake.output.gene_log2_txt}.md5
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

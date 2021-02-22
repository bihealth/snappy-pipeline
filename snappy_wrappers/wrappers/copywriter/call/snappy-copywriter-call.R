library(magrittr)

# WARING- the order of library loading is important!
library( matrixStats )
library( Rsamtools )
library( DNAcopy )
library( CopywriteR )
library( CGHcall )

postProcess <- function(destination.folder="./work/donor",
                        donorID="donor",
                        fullID="pk-donor",
                        cghcall.params=list(
                            inter=c( -0.4, 0.4 ),
                            nclass=5, divide=5,
                            robustsig="yes",
                            cellularity=0.6,
                            ncpus=4),
                        genomeRelease="hg19",
                        features=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75,
                        plot.params=list( 
                            genes=NULL,
                            main="",
                            ylim=c( -3, 4 ) ),
                        mapper="bwa") {
    if( is.null( destination.folder ) ||
        !is.character( destination.folder ) ||
        length( destination.folder ) != 1 ||
        any( is.na( destination.folder ) ) ) {
        stop( "Illegal destination folder" )
    }
    if( !dir.exists( destination.folder ) ) {
        dir.create( destination.folder, showWarnings=FALSE, recursive=TRUE, mode="0750" )
        if( !dir.exists( destination.folder ) ) {
            stop( "Can't create destination folder" )
        }
    }
    
    if( is.null( cghcall.params ) || !is.list( cghcall.params ) ||
        is.null( names( cghcall.params ) ) ||
        !all( c( "inter", "nclass", "divide", "robustsig", "cellularity", "ncpus" ) %in% names( cghcall.params ) ) ) {
        stop( "Missing parameter(s) for CGH calls" )
    }
    
    if( is.null( plot.params ) || !is.list( plot.params ) ||
        is.null( names( plot.params ) ) ||
        !all( c( "genes", "main", "ylim" ) %in% names( plot.params ) ) ) {
        stop( "Missing parameter(s) for plotting" )
    }
    
    # load the input parameters to have the prefix value on hand
    load(file.path(destination.folder, "CNAprofiles/input.Rdata"))
    sampleName <- paste(paste("log2", make.names(basename(inputStructure$sample.control$samples)), sep="."), collapse=".vs.")
    
    load( file.path( destination.folder, "CNAprofiles/segment.Rdata" ) )
    segment.CNA.object <- subset.CNA.object( segment.CNA.object, sampleName )
    if( length.CNA.object( segment.CNA.object ) != 1 ) {
        stop( "Unexpected number of contrasts" )
    }

    # Avoid using scientific notation for locations
    options( scipen=999 )

    # Prepare segments & bin files for output
    segData <- getSegData( segment.CNA.object, allowedChrom=inputStructure$chromosomes )
    assertthat::assert_that( length( unique( segData$ID ) )==1 )
    write.table( segData %>% dplyr::mutate( ID=donorID ) %>% dplyr::select( -startRow, -endRow ),
                 file.path( destination.folder, "out",
                            paste( mapper, ".copywriter.", fullID, "_segments.txt", sep="" ) ),
                 sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )

    binData <- getBinData( segment.CNA.object, binSize=inputStructure$bin.size, allowedChrom=inputStructure$chromosomes )
    assertthat::assert_that( ncol( binData )==6 )
    write.table( binData %>% dplyr::select( -iBin ),
                 file.path( destination.folder, "out",
                            paste( mapper, ".copywriter.", fullID, "_bins.txt", sep="" ) ),
                 sep="\t", row.names=FALSE, quote=FALSE,
                 col.names=c( "chrom", "maploc", donorID, "loc.start", "loc.end" ) )

    # Plots
    plot.dir <- file.path( destination.folder, "plots" ) 
    dir.create( plot.dir, showWarnings=FALSE, recursive=FALSE, mode="0750" )
    chrLen <- getChromosomeLengths( inputStructure$sample.control$samples[1], inputStructure$chromosomes )
    plotOneSample( segment.CNA.object, chrLen, binSize=inputStructure$bin.size,
                   sample=sampleName,
                   main=donorID,
                   genes=plot.params$genes,
                   plot.dir=plot.dir,
                   ylim=plot.params$ylim )

    # Create the calls using CGHcall
    seg <- to_cghSeg( segment.CNA.object, allowedChrom=inputStructure$chromosomes, binSize=inputStructure$bin.size, genomeRelease=genomeRelease )
    sgn <- postsegnormalize( seg, inter=cghcall.params$inter )
    cll <- CGHcall( sgn, nclass=cghcall.params$nclass, robustsig=cghcall.params$robustsig, 
                    cellularity=cghcall.params$cellularity, ncpus=cghcall.params$ncpus )
    res <- ExpandCGHcall( cll, sgn, divide=cghcall.params$divide, memeff=FALSE )

    # to_cghSeg set chromosome names to 1:24 for compatibility with CGHcall
    # the actual chromosome names are restored in res
    tmp <- pData(featureData(res))
    tmp$Chromosome <- chrLen$Chromosome[tmp$Chromosome]
    pData(featureData(res)) <- tmp
    
    save( res, file=file.path( destination.folder, "segment.RData" ) )

    # Output the calls
    tmp <- geneCalls( segData, res, binSize=inputStructure$bin.size, type="call", tx_obj=features )
    i   <- which( colnames( tmp ) == sampleName )
    colnames( tmp )[i] <- donorID
    write.table( tmp,
                 file=file.path( destination.folder, "out",
                                 paste( mapper, ".copywriter.", fullID, "_gene_call.txt", sep="" ) ),
                 sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
    tmp <- geneCalls( segData, res, binSize=inputStructure$bin.size, type="log2", tx_obj=features )
    i   <- which( colnames( tmp ) == sampleName )
    colnames( tmp )[i] <- donorID
    write.table( tmp,
                 file=file.path( destination.folder, "out",
                                 paste( mapper, ".copywriter.", fullID, "_gene_log2.txt", sep="" ) ),
                 sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE )
}

subset.CNA.object <- function( segment.CNA.object, contrasts ) {
    if( missing( segment.CNA.object ) || is.null( segment.CNA.object ) ||
        class( segment.CNA.object ) != "DNAcopy" ) {
        stop( "Illegal CNA object" )
    }
    if( missing( contrasts ) || is.null( contrasts ) ||
        !any( contrasts %in% as.character( segment.CNA.object$output$ID ) ) ) {
        tmp <- list( data=CNA( matrix( 0, nrow=nrow( segment.CNA.object$data ), ncol=0 ),
                               segment.CNA.object$data$chrom,
                               segment.CNA.object$data$maploc,
                               attr( segment.CNA.object$data, "data.type" ),
                               sampleid=NULL, presorted=FALSE ),
                     output=data.frame( ID="", chrom=0, loc.start=0, len.end=0, num.mark=0, seg.mean=0 )[-1,],
                     segRows=data.frame( startRow=0, endRow=0 )[-1,],
                     call=segment.CNA.object$call )
        class( tmp ) <- "DNAcopy"
        return( tmp )
    }
    contrasts <- contrasts[contrasts %in% colnames( segment.CNA.object$data )]
    i   <- which( as.character( segment.CNA.object$output$ID ) %in% contrasts )
    tmp <- list( data=CNA( as.matrix( segment.CNA.object$data )[,contrasts],
                           segment.CNA.object$data$chrom,
                           segment.CNA.object$data$maploc,
                           attr( segment.CNA.object$data, "data.type" ),
                           sampleid=contrasts, presorted=FALSE ),
                 output=segment.CNA.object$output[i,],
                 segRows=segment.CNA.object$segRows[i,],
                 call=segment.CNA.object$call )
    class( tmp ) <- "DNAcopy"
    tmp
}

length.CNA.object <- function( segment.CNA.object ) {
    sum( !(colnames( segment.CNA.object$data ) %in% c( "chrom", "maploc" )) )
}

getSegData <- function( segment.CNA.object, allowedChrom ) {
    segData <- cbind( segment.CNA.object[["output"]], segment.CNA.object$segRows ) %>%
        dplyr::mutate( loc.start=ceiling( loc.start ), loc.end=floor( loc.end ) )
    
    segData$chrom <- factor(allowedChrom[segData$chrom], levels=allowedChrom)
    
    segData
}

getBinData <- function( segment.CNA.object, binSize=20000, allowedChrom ) {
    binData <- as.data.frame( segment.CNA.object[["data"]] ) %>%
        dplyr::mutate( loc.start=(2*maploc - (binSize-1))/2,
                       loc.end=(2*maploc + (binSize-1))/2 )
    binData$iBin <- 1:nrow( binData )
    assertthat::assert_that( all( abs( binData$loc.start - floor( binData$loc.start ) ) < 0.001 ) )
    assertthat::assert_that( all( abs( binData$loc.end   - floor( binData$loc.end   ) ) < 0.001 ) )
    
    binData$chrom <- factor(allowedChrom[binData$chrom], levels=allowedChrom)

    binData
}

to_cghSeg <- function( segment.CNA.object, binSize=20000, allowedChrom, genomeRelease="hg19" ) {
    binData <- getBinData( segment.CNA.object, binSize=binSize, allowedChrom=allowedChrom )
    binData <- binData %>%
        dplyr::mutate( Feature=sprintf( "%s:%.0f-%.0f", chrom, loc.start, loc.end ) )
    valCols <- which( !(colnames( binData ) %in% c( "chrom", "maploc", "loc.start", "loc.end", "iBin", "Feature" )) )

    segData <- getSegData( segment.CNA.object, allowedChrom=allowedChrom )

    cghData <- matrix( NA, nrow=nrow( binData ), ncol=length( valCols ),
                       dimnames=list( binData$Feature, colnames( binData )[valCols] ) )

    for( i in 1:nrow( segData ) ) {
        j <- segData$startRow[i]:segData$endRow[i]
        cghData[j,segData$ID[i]] <- segData$seg.mean[i]
    }

    i <- which( rowSums( is.na( binData[,valCols,drop=FALSE] ) ) == 0 )

    sampleNames <- colnames( binData )[valCols]
    df.features <- AnnotatedDataFrame( data.frame( Chromosome=as.numeric(binData$chrom[i]),
                                                   Start=binData$loc.start[i], End=binData$loc.end[i],
                                                   row.names=binData$Feature[i],
                                                   check.names=FALSE ),
                                       varMetadata=data.frame( labelDescription=c( "Chromosomal position",
                                                                                   "Basepair position start",
                                                                                   "Basepair position end" ),
                                                               row.names=c( "Chromosome", "Start", "End" ) ),
                                       dimLabels=c( "featureNames", "colnames" ) )
    df.pheno    <- AnnotatedDataFrame( data.frame( x=rep( 0, length( sampleNames ) ),
                                                   row.names=sampleNames )[,-1],
                                       varMetadata=as.data.frame( matrix( NA, nrow=0, ncol=1, 
                                                                  dimnames=list( NULL, "labelDescription" ) ) ) )

    rownames( binData ) <- binData$Feature
    binData <- as.matrix( binData[rownames( cghData ),colnames( cghData ),drop=FALSE] )

    new( "cghSeg", copynumber=binData[i,,drop=FALSE], segmented=cghData[i,,drop=FALSE],
         featureData=df.features, phenoData=df.pheno, annotation=genomeRelease )
}

to_GRanges <- function( segData, binSize=20000, 
                        tx_obj=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) {
    seg <- GRanges( seqnames=Rle( segData$chrom ),
                    ranges=IRanges( start=segData$loc.start, end=segData$loc.end ),
                    seqinfo=GenomeInfoDb::seqinfo(tx_obj) )

    mcols( seg ) <- segData[,c( 1, 5:ncol( segData ) )]
    seg
}

segmentCalls <- function( segData, res, binSize=20000 ) {
    tmp <- calls( res )
    tmp <- data.frame( tmp,
                       chrom=sub( ":.*", "", rownames( tmp ) ),
                       loc.start=as.integer( sub( ".*:", "", sub( "-.*", "", rownames( tmp ) ) ) ),
                       loc.end=as.integer( sub( ".*-", "", rownames( tmp ) ) ) )

    segData$call <- NA
    for( i in 1:nrow( segData ) ) {
        j1 <- which( tmp$chrom==segData$chrom[i] &
                     tmp$loc.start<=segData$loc.start[i] & segData$loc.start[i]<tmp$loc.end )
        j2 <- which( tmp$chrom==segData$chrom[i] &
                     tmp$loc.start<=segData$loc.end[i] & segData$loc.end[i]<tmp$loc.end )
        if( length( j1 )==1 && length( j2 )==1 ) {
            v <- unique( tmp[j1:j2,segData$ID[i]] )
            if( length( v ) == 1 ) segData$call[i] <- v
        }
    }

    segData
}

geneCalls <- function( segData, res, binSize=20000, type=c( "call", "log2" ),
                       ann=org.Hs.eg.db::org.Hs.eg.db, 
                       tx_obj=EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75) {
    type <- match.arg( type )

    geneNames <- AnnotationDbi::select( ann, columns="SYMBOL",
                                        keys=AnnotationDbi::keys( ann, keytype="ENTREZID" ), keytype="ENTREZID" )
    i <- which( is.na( geneNames$ENTREZID ) | is.na( geneNames$SYMBOL ) |
                geneNames$ENTREZID=="" | as.numeric( geneNames$ENTREZID ) <= 0 |
                geneNames$SYMBOL=="" )
    if( length( i ) > 0 ) geneNames <- geneNames[-i,,drop=FALSE]
    if( nrow( geneNames ) == 0 )
        stop( "No genes" )

    genes <- GenomicFeatures::genes( tx_obj )

    # Retain only 1-1 between ENSEMBL & ENTREZ ids, & set gene_id to entrezid
    if (is(tx_obj, "EnsDb")) {
        i <- sapply(genes$entrezid, function(l) length(l)==1 && !is.na(l))
        genes <- genes[i,,drop=FALSE]
        tmp <- GenomicRanges::mcols(genes)
        tmp$gene_id <- as.character(unlist(tmp$entrezid))
        GenomicRanges::mcols(genes) <- tmp
    }

    i <- which( mcols( genes )$gene_id %in% geneNames$ENTREZID )
    if( length( i ) == 0 )
        stop( "No common genes" )
    genes <- genes[i,,drop=FALSE]

    seg <- to_GRanges( segmentCalls( segData, res, binSize=binSize ), binSize=binSize, tx_obj=tx_obj )

    sampleNames <- unique( segData$ID )
    rslt <- matrix( NA, nrow=length( genes ), ncol=length( sampleNames ),
                    dimnames=list( mcols( genes )$gene_id, sampleNames ) )
    for( sample in sampleNames ) {
        s <- seg[mcols( seg )$ID==sample,]
        i <- as.matrix( findOverlaps( genes, s ) )
        i <- split( i[,"subjectHits"], i[,"queryHits"] )
        i <- unlist( i[lengths( i )==1] )
        if( type == "call" )
            rslt[as.numeric( names( i ) ),sample] <- mcols( s )$call[i]
        if( type == "log2" )
            rslt[as.numeric( names( i ) ),sample] <- mcols( s )$seg.mean[i]
    }

    i <- match( mcols( genes )$gene_id, geneNames$ENTREZID )
    geneNames <- geneNames[i,]
    data.frame( Hugo_Symbol=geneNames$SYMBOL,
                Entrez_Gene_Id=geneNames$ENTREZID,
                rslt )
}

getChromosomeLengths <- function( bam, allowedChroms ) {
    chrLen <- scanBamHeader( bam )[[1]]$targets
    chrLen <- data.frame( Chromosome=names( chrLen ),
                       Length=chrLen,
                         row.names=NULL, stringsAsFactors=FALSE )
    chrLen <- chrLen[chrLen$Chromosome %in% allowedChroms,]
    chrLen$Chromosome <- factor( chrLen$Chromosome, levels=allowedChroms )
    chrLen[order( chrLen$Chromosome ),]
}

plotOneSample <- function( segment.CNA.object, chrLen, binSize=20000,
                           sample="log2.tumor.bam.vs.log2.normal.bam",
                           genes=NULL, plot.dir=".",
                           main="donor", ylim=c( -3, 4 ) ) {
    segData <- getSegData( segment.CNA.object, levels(chrLen$Chromosome) )
    segData <- segData %>%
        dplyr::filter( ID==sample )

    binData <- getBinData( segment.CNA.object, binSize=binSize, levels(chrLen$Chromosome) )
    binData <- binData[,c( "chrom", "maploc", sample )]
    colnames( binData )[3] <- "value"

    for( i in levels(chrLen$Chromosome)) {
        if( capabilities( "png" ) ) {
            png( file.path( plot.dir, paste( i, "png", sep="." ) ),
                 height=720, width=1440 )
        } else {
            bitmap( file.path( plot.dir, paste( i, "png", sep="." ) ),
                    height=720, width=1440, unit="px", type="png16m" )
        }
        plotOneChrom( binData, segData, chrLen, myChrom=i, genes=genes, ylim=ylim )
                      title( main=paste( main, i, sep=" - Chromosome " ) )
        dev.off()
    }

    if( capabilities( "png" ) ) {
        png( file.path( plot.dir, "All.png" ),
                        height=720, width=1440 )
    } else {
        bitmap( file.path( plot.dir, "All.png" ),
                height=720, width=1440, unit="px", type="png16m" )
    }
    plotAllChrom( binData, segData, chrLen, genes=genes, ylim=ylim )
    title( main=main )
    dev.off()
}

plotOneChrom <- function( binData, segData, chrLen,
                          myChrom="1", genes=NULL,
                          ylim=c( -3, 4 ) ) {
    ylim <- sort( ylim )

    chrLen <- chrLen$Length[chrLen$Chromosome==myChrom]

    segData <- segData %>%
        dplyr::filter( chrom==myChrom ) %>%
        dplyr::mutate( center=0.5*(loc.start + loc.end) ) %>%
        dplyr::mutate( cropped=seg.mean ) %>%
        dplyr::mutate( cropped=replace( cropped, seg.mean < ylim[1], ylim[1] ) ) %>%
        dplyr::mutate( cropped=replace( cropped, seg.mean > ylim[2], ylim[2] ) ) %>%
        dplyr::mutate( col=2 ) %>%
        dplyr::mutate( col=replace( col, seg.mean!=cropped, 5 ) ) %>%
        dplyr::arrange( center )

    binData <- binData %>%
        dplyr::filter( chrom==myChrom ) %>%
        dplyr::filter( !is.na( value ) ) %>%
        dplyr::mutate( cropped=value ) %>%
        dplyr::mutate( cropped=replace( cropped, value < ylim[1], ylim[1] ) ) %>%
        dplyr::mutate( cropped=replace( cropped, value > ylim[2], ylim[2] ) ) %>%
        dplyr::mutate( pch="." ) %>%
        dplyr::mutate( pch=replace( pch, value!=cropped, "+" ) )

    plot( c( 1, chrLen ), ylim+c( 0,!is.null( genes ) ), type="n",
          xlab="Position", ylab="CNV (log2 scale)" )

    s <- madDiff( binData$value )
    abline( h=3*s*c( -1, 1 ), col=1, lty=3 )
    text( chrLen, ylim[2]+(!is.null( genes )), 
          labels=sprintf( "MAD: %6.3f", s ), adj=c( 1, 1 ) )

    points( binData$maploc, binData$cropped, 
           col=1, pch=binData$pch, cex=1.5 )

    for( i in 1:nrow( segData ) ) {
        lines( c( segData$loc.start[i], segData$loc.end[i] ),
               rep( segData$cropped[i], 2 ), col=segData$col[i], lwd=2 )
    }

    if( !is.null( genes ) ) {
        if( sum( genes$chromosome_name==myChrom ) > 0 ) {
            genes <- genes[genes$chromosome_name==myChrom,]
            genes$center <- 0.5*(genes$start_position + genes$end_position)
            i            <- findInterval( genes$center, segData$loc.start )
            genes$ypos   <- segData$cropped[i]
                        
            points( genes$center, genes$ypos, pch=16, col=3 )
            text( genes$center, ylim[2]+0.8, labels=genes$hgnc_symbol,
                  adj=c( 0, -0.5 ), srt=90 )
        }
    }
}

plotAllChrom <- function( binData, segData, chrLen,
                          genes=NULL,
                          ylim=c( -3, 4 ) ) {
    ylim <- sort( ylim )

    chrLen$CumSum <- c( 0, cumsum( as.numeric( chrLen$Length[-nrow( chrLen )] ) ) )
    chrLen$Center <- chrLen$CumSum + chrLen$Length/2

    segData <- segData %>%
        dplyr::left_join( chrLen, by=c( "chrom"="Chromosome" ) ) %>%
        dplyr::mutate( loc.start=loc.start + CumSum,
                       loc.end=loc.end + CumSum ) %>%
        dplyr::mutate( center=0.5*(loc.start + loc.end) ) %>%
        dplyr::mutate( cropped=seg.mean ) %>%
        dplyr::mutate( cropped=replace( cropped, seg.mean < ylim[1], ylim[1] ) ) %>%
        dplyr::mutate( cropped=replace( cropped, seg.mean > ylim[2], ylim[2] ) ) %>%
        dplyr::mutate( col=2 ) %>%
        dplyr::mutate( col=replace( col, seg.mean!=cropped, 5 ) ) %>%
        dplyr::arrange( center )

    binData <- binData %>%
        dplyr::left_join( chrLen, by=c( "chrom"="Chromosome" ) ) %>%
        dplyr::mutate( maploc=maploc + CumSum ) %>%
        dplyr::filter( !is.na( value ) ) %>%
        dplyr::mutate( cropped=value ) %>%
        dplyr::mutate( cropped=replace( cropped, value < ylim[1], ylim[1] ) ) %>%
        dplyr::mutate( cropped=replace( cropped, value > ylim[2], ylim[2] ) ) %>%
        dplyr::mutate( pch="." ) %>%
        dplyr::mutate( pch=replace( pch, value!=cropped, "+" ) )

    plot( c( 0, max( chrLen$CumSum ) ), ylim+c( 0, !is.null( genes ) ), type="n",
          xaxt="n", xlab="Position", ylab="CNV (log2 scale)" )
    abline( v=chrLen$CumSum[-1], col=1, lty=3 )
    axis( 1, at=chrLen$Center, label=as.character( chrLen$Chromosome ) )

    s <- madDiff( binData$value )
    abline( h=3*s*c( -1, 1 ), col=1, lty=3 )
    text( max( chrLen$CumSum ), ylim[2]+(!is.null( genes )), 
          labels=sprintf( "MAD: %6.3f", s ), adj=c( 1, 1 ) )

    points( binData$maploc, binData$cropped, 
            col=1, pch=binData$pch, cex=1.5 )

    for( i in 1:nrow( segData ) ) {
        lines( c( segData$loc.start[i], segData$loc.end[i] ),
               rep( segData$cropped[i], 2 ), col=segData$col[i], lwd=2 )
    }

    if( !is.null( genes ) ) {
        i             <- match( as.character( genes$chromosome_name ),
                                as.character( chrLen$Chromosome ) )
        genes$start_position <- genes$start_position + chrLen$CumSum[i]
        genes$end_position   <- genes$end_position   + chrLen$CumSum[i]
        genes$center         <- 0.5*(genes$start_position + genes$end_position)
        i             <- findInterval( genes$center, segData$loc.start )
        genes$ypos    <- segData$cropped[i]

        points( genes$center, genes$ypos, pch=16, col=3 )
    }
}

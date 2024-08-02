library(magrittr)

# WARING- the order of library loading is important!
library( matrixStats )
library( Rsamtools )
library( CopywriteR )

runCopywriteR <- function( normal.bam, tumor.bam,
                    destination.folder="./work/donor",
                    binSize=20000,
                    donorID="donor",
                    fullID="pk-donor",
                    copywriter.params=list(
                        reference.folder="work/hg19_20kb",
                        capture.regions.file="/fast/projects/cubit/current/static_data/exome_panel/Agilent/SureSelect_Human_All_Exon_V5/GRCh37/Regions.bed",
                        workers=8 )
                    ) {
    if( missing( normal.bam ) || missing( tumor.bam ) ||
        is.null( normal.bam ) || is.null( tumor.bam ) ||
        !is.character( normal.bam ) || !is.character( tumor.bam ) ||
        length( normal.bam ) != 1 || length( tumor.bam ) != 1 ||
        any( is.na( normal.bam ) ) || any( is.na( tumor.bam ) ) ) {
        stop( "Illegal normal or tumor file names" )
    }
    if( !file.exists( normal.bam ) || !file.exists( tumor.bam ) ) {
        stop( "Missing normal or tumor file" )
    }

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
    
    if( is.null( copywriter.params$reference.folder ) ||
        !is.character( copywriter.params$reference.folder ) ||
        length( copywriter.params$reference.folder ) != 1 ||
        any( is.na( copywriter.params$reference.folder ) ) ) {
        stop( "Illegal reference folder" )
    }
    if( !file.exists( file.path( copywriter.params$reference.folder, "GC_mappability.rda" ) ) ||
        !file.exists( file.path( copywriter.params$reference.folder, "blacklist.rda" ) ) ) {
        stop( "Missing GC, mappability or blacklist files" )
    }

    if( is.null( copywriter.params$capture.regions.file ) ||
        !is.character( copywriter.params$capture.regions.file ) ||
        length( copywriter.params$capture.regions.file ) != 1 ||
        any( is.na( copywriter.params$capture.regions.file ) ) ) {
        stop( "Illegal capture regions filefolder" )
    }
    if( !file.exists( copywriter.params$capture.regions.file ) ) {
        stop( "Missing capture regions file" )
    }

    if( is.null( binSize ) ||
        !is.numeric( binSize ) ||
        length( binSize ) != 1 ||
        any( is.na( binSize ) || is.infinite( binSize ) ) ) {
        stop( "Illegal bin size length" )
    }
    if( binSize <= 0 ) {
        stop( "Bin size must be strictly positive" )
    }

    if( is.null( copywriter.params ) || !is.list( copywriter.params ) ||
        is.null( names( copywriter.params ) ) ||
        !all( c( "reference.folder", "capture.regions.file", "workers" ) %in% names( copywriter.params ) ) ) {
        stop( "Missing parameter(s) for copywriter" )
    }

    # Workaround lack of png capabilities on the cluster (use bitmap)
    assignInNamespace( ".tng", value=my.tng, ns="CopywriteR" )

    bp.param <- SnowParam( workers=copywriter.params$workers, type="SOCK" )

    dir.create( destination.folder, showWarnings=FALSE, recursive=TRUE, mode="0750" )
    if( file.exists( file.path( destination.folder, "CNAprofiles" ) ) ) {
        unlink( file.path( destination.folder, "CNAprofiles" ), recursive=TRUE )
    }

    set.seed( 1234567 )
    CopywriteR( sample.control=cbind( c( tumor.bam, normal.bam ),
                                      c( normal.bam, normal.bam ) ),
                destination.folder=destination.folder,
                reference.folder=copywriter.params$reference.folder,
                capture.regions.file=copywriter.params$capture.regions.file,
                bp.param=bp.param )

    plotCNA( destination.folder=destination.folder )
}

my.tng <- function(df, use, correctmappa = TRUE, plot = NULL, verbose = TRUE) {
    #tests
    if (!is.logical(use) && length(use) == nrow(df)) 
        stop("use should be logicval vector with same size as df")
    #df colums?
    
    if (!is.null(plot)) {
        if (!is.logical(plot)) {
            if (verbose) 
                cat("Plotting to file", plot, "\n")
            if( capabilities( "png" ) ) {
                png(plot, width = 700, height = 1400)
            } else {
                bitmap( plot, type="png16m", width=700, height=1400, unit="px" )
            }
            par(mfrow = c(2, 1))
            on.exit(dev.off())
            plot <- TRUE
        } else if (plot) {
            par(mfrow = c(2, 1))
        }
    }

    #exclude contains the points to exclude in the 
    #fitting (usually sex chromosomes and blacklisted regions)
    #gc fits also excludes the low mappability data

    #correct gc using double lowess
    gcuse <- (use & !is.na(df$mappa) & df$mappa > 0.8 & !is.na(df$gc) & df$gc > 0)
    rough <- loess(count ~ gc, data = df, subset = gcuse, span = 0.03)
    i <- seq(0, 1, by = 0.001)
    final <- loess(predict(rough, i) ~ i, span = 0.3)
    normv <- predict(final, df$gc)
    df$countgcloess <- df$count/(normv/median(normv, na.rm = TRUE))

    if (plot) {
        plot(count ~ gc, data = df, subset = gcuse,
             ylim = quantile(df$count[gcuse], c(1e-04, 0.999)), xlim = c(0, 1),
             pch = ".")
        points(count ~ gc, data = df, subset = !gcuse, col = rgb(1, 0, 0, 0.3),
               pch = ".")
        lines(i, predict(rough, i), col = "green")
        points(df$gc, normv, col = "red", pch = ".")
    }

    #correct mappa using linear function that intercepts zero
    #if(correctmappa) {
    #mappause <- (use & !is.na(df$mappa))
    #lm(countgcloess~0+mappa, data=df, subset=mappause) ->fll
    #if(verbose) print(summary(fll))

    #if (plot) {
    #  plot(countgcloess ~ mappa, data=df, subset=mappause,
    #       ylim=quantile(df$countgcloess, c(0.0001, .999), na.rm=T), pch=".")
    #  points(countgcloess ~ mappa, data=df, subset=!mappause,
    #         col=rgb(1,0,0,.3), pch=".")
    #  abline(0, fll$coef, col=2)
    #}

    #correct mappa using double lowess -> paired end sequencing
    if (correctmappa) {
        mappause <- (use & !is.na(df$mappa))
        rough <- loess(countgcloess ~ mappa, data = df, subset = mappause,
                       span = 0.03)
        i <- seq(0, 1, by = 0.001)
        final <- loess(predict(rough, i) ~ i, span = 0.3)
        normv <- predict(final, df$mappa)
        df$countgcmappaloess <- df$countgcloess/(normv/median(normv,
                                                              na.rm = TRUE))

        if (plot) {
            plot(countgcloess ~ mappa, data = df, subset = mappause,
                 ylim = quantile(df$countgcloess[mappause], c(1e-04, 0.999),
                                 na.rm = TRUE), xlim = c(0, 1), pch = ".")
            points(countgcloess ~ mappa, data = df, subset = !mappause,
                   col = rgb(1, 0, 0, 0.3), pch = ".")
            lines(i, predict(rough, i), col = "green")
            points(df$mappa, normv,
                   subset = (!is.na(df$mappa) && !is.na(normv)),
                   col = "red", pch = ".")
        }

        return(log2(df$countgcmappaloess/median(df$countgcmappaloess[use],
               na.rm = TRUE)))
    } else {
        #corerct agains median value (exluding sex chr)
        log2(df$countgcloess/median(df$countgcloess[use], na.rm = TRUE))
    }
}


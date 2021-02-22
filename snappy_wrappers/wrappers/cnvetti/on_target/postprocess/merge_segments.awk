BEGIN {
    FS = "\t";
    print "ID", "chrom", "loc.start", "loc.end", "num.mark", "seg.mean";
}
NR==1 {CHROM=$3 ; SEGSTART=$4 ; SEGEND=$5 ; SEGMEAN=$6; SEGCOUNTER=0}
{
    if ($6 == SEGMEAN)
    {
        SEGEND=$5;
        SEGCOUNTER++;
    }
    else
    {
        print $1, CHROM, SEGSTART, SEGEND, SEGCOUNTER, SEGMEAN;
        CHROM=$3 ; SEGSTART=$4 ; SEGEND=$5 ; SEGMEAN=$6; SEGCOUNTER=1;
    }
}

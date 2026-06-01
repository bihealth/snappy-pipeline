# -*- coding: utf-8 -*-
"""Wrapper for actually running ASCAT"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
# -------------------------------------------------------------------------------------------------
# Build input files
#
grep -v '^.chrs' {snakemake.input.baf_tumor} > $TMPDIR/baf_tumor.raw.txt
grep -v '^.chrs' {snakemake.input.cnv_tumor} > $TMPDIR/cnv_tumor.raw.txt
grep -v '^.chrs' {snakemake.input.baf_normal} > $TMPDIR/baf_normal.raw.txt
grep -v '^.chrs' {snakemake.input.cnv_normal} > $TMPDIR/cnv_normal.raw.txt

comm -12 \
    <(  \
        comm -12 \
            <(cut -f 1 $TMPDIR/baf_tumor.raw.txt | sort) \
            <(cut -f 1 $TMPDIR/cnv_tumor.raw.txt | sort) \
    ) \
    <( \
        comm -12 \
            <(cut -f 1 $TMPDIR/baf_normal.raw.txt | sort) \
            <(cut -f 1 $TMPDIR/cnv_normal.raw.txt | sort) \
    ) \
| grep -v '^$' \
> $TMPDIR/ids.txt

for name in baf_tumor baf_normal cnv_tumor cnv_normal; do
    echo -e "\tchrs\tpos\t{args[tumor_library_name]}" > $TMPDIR/$name.txt
    grep -w -F -f $TMPDIR/ids.txt $TMPDIR/$name.raw.txt >> $TMPDIR/$name.txt
done

# -------------------------------------------------------------------------------------------------
# Write ASCAT script and call R
#
cat <<"EOF" > $TMPDIR/run_ascat.R
library(ASCAT)

GC_CONTENT_FILE = "/fast/users/mholtgr/Data/ASCAT/GC_AffySNP6_102015.txt";

ascat.bc = ascat.loadData(
    "cnv_tumor.txt",
    "baf_tumor.txt",
    "cnv_normal.txt",
    "baf_normal.txt");

#ascat.bc = ascat.GCcorrect(ascat.bc, GC_CONTENT_FILE)
ascat.plotRawData(ascat.bc)
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

# Write out results.
write.table(ascat.output$nA, "{args[tumor_library_name]}_na.txt");
write.table(ascat.output$nB, "{args[tumor_library_name]}_nb.txt");
write.table(ascat.output$goodnessOfFit, "{args[tumor_library_name]}_goodness_of_fit.txt");
write.table(ascat.output$ploidy, "{args[tumor_library_name]}_ploidy.txt");
write.table(ascat.output$psi, "{args[tumor_library_name]}_psi.txt");
write.table(ascat.output$segments, "{args[tumor_library_name]}_segments.txt");
write.table(ascat.output$segments_raw, "{args[tumor_library_name]}_segments_raw.txt");
write.table(ascat.output$aberrantcellfraction, "{args[tumor_library_name]}_aberrantcellfraction.txt");
EOF

out_dir=$(readlink -f $(dirname {snakemake.output.done}))

cd $TMPDIR
R --vanilla < run_ascat.R

# -------------------------------------------------------------------------------------------------
# Move out output
#

for suffix in .ASCATprofile.png .ASPCF.png .germline.png .rawprofile.png .sunrise.png .tumour.png \
              _na.txt _nb.txt _goodness_of_fit.txt _ploidy.txt _segments.txt _segments_raw.txt \
              _psi.txt _aberrantcellfraction.txt; do
    if [ -f "{args[tumor_library_name]}$suffix" ]; then
        cp \
            {args[tumor_library_name]}$suffix \
            $out_dir
    else
        echo "WARNING {args[tumor_library_name]}$suffix -- File does not exist"
    fi

done
"""
)

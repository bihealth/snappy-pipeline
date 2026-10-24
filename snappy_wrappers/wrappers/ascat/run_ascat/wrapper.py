# -*- coding: utf-8 -*-
"""Wrapper for actually running ASCAT"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# ASCAT names its files after the tumor library; the declared outputs are the only ones kept.
# The R script runs in $TMPDIR, so the targets are absolute.
ascat_files = {
    "goodness_of_fit": "goodness_of_fit.txt",
    "ploidy": "ploidy.txt",
    "purity": "aberrantcellfraction.txt",
    "segments": "segments.txt",
    "segments_raw": "segments_raw.txt",
    "plot_tumour": "tumour.png",
    "plot_germline": "germline.png",
    "plot_aspcf": "ASPCF.png",
}
copy_results = "\n".join(
    f"cp {args['tumor_library']}.{suffix} {os.path.abspath(snakemake.output[key])}"
    for key, suffix in ascat_files.items()
)

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
    <(
        comm -12 \
            <(cut -f 1 $TMPDIR/baf_tumor.raw.txt | sort) \
            <(cut -f 1 $TMPDIR/cnv_tumor.raw.txt | sort)
    ) \
    <(
        comm -12 \
            <(cut -f 1 $TMPDIR/baf_normal.raw.txt | sort) \
            <(cut -f 1 $TMPDIR/cnv_normal.raw.txt | sort)
    ) \
    | grep -v '^$' \
        > $TMPDIR/ids.txt

for name in baf_tumor baf_normal cnv_tumor cnv_normal; do
    echo -e "\tchrs\tpos\t{args[tumor_library]}" > $TMPDIR/$name.txt
    grep -w -F -f $TMPDIR/ids.txt $TMPDIR/$name.raw.txt >> $TMPDIR/$name.txt
done

# -------------------------------------------------------------------------------------------------
# Write ASCAT script and call R
#
cat << "EOF" > $TMPDIR/run_ascat.R
library(ASCAT)

ascat.bc = ascat.loadData(
    "cnv_tumor.txt",
    "baf_tumor.txt",
    "cnv_normal.txt",
    "baf_normal.txt");

ascat.plotRawData(ascat.bc)
ascat.plotRawData(ascat.bc)
ascat.bc = ascat.aspcf(ascat.bc)
ascat.plotSegmentedData(ascat.bc)
ascat.output = ascat.runAscat(ascat.bc)

# Write out results.
write.table(ascat.output$goodnessOfFit, "{args[tumor_library]}.goodness_of_fit.txt");
write.table(ascat.output$ploidy, "{args[tumor_library]}.ploidy.txt");
write.table(ascat.output$segments, "{args[tumor_library]}.segments.txt");
write.table(ascat.output$segments_raw, "{args[tumor_library]}.segments_raw.txt");
write.table(ascat.output$aberrantcellfraction, "{args[tumor_library]}.aberrantcellfraction.txt");
EOF

cd $TMPDIR
R --vanilla < run_ascat.R

# -------------------------------------------------------------------------------------------------
# Keep the declared output
#
{copy_results}
"""
)

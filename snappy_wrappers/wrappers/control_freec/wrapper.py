# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Control-FreeC: Snakemake wrapper.py"""

from snakemake import shell

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

window_str = ""
w = args.get("window_size")
if w >= 0:
    window_str = "window = {}".format(w)

mappability = snakemake.input.get("mappability", "")

contamination_adjustment = "TRUE" if args["contamination_adjustment"] else "FALSE"
unique_match = "TRUE" if args["unique_match"] else "FALSE"

shell.executable("/bin/bash")

ShellWrapper(snakemake).run(
    r"""
set -x

# write out configuration -------------------------------------------------------------------------
export output_dir=$TMPDIR/freec_out

mkdir -p $TMPDIR/tmp.d $output_dir

cat << EOF > $TMPDIR/tmp.d/freec.ini
# http://boevalab.com/FREEC/tutorial.html
[general]

maxThreads = {snakemake.threads}
outputDir = $output_dir

## path to sambamba (faster BAM file reading)
sambamba = sambamba

chrLenFile = {snakemake.input.chrlenfile}
ploidy = {args[ploidy]}

breakPointThreshold = {args[break_point_threshold]}
coefficientOfVariation = {args[coefficient_of_variation]}
{window_str}

# normal cell contamination of the tumor sample
contaminationAdjustment = {contamination_adjustment}
contamination = {args[contamination]}
minimalSubclonePresence = {args[minimal_subclone_presence]}

numberOfProcesses = {snakemake.threads}

$(if [[ "{args[path_mappability_enabled]}" == True ]]; then
    echo gemMappabilityFile = {mappability}
fi)

uniqueMatch = {unique_match}
forceGCcontentNormalization = {args[force_gc_content_normalization]}

minCNAlength = {args[min_cna_length]}

BedGraphOutput=TRUE

# this makes CNV predictions on Gonosomes unreliable
sex={args[sex]}

# this sets the degree of polynomial; currently the step of estimating it seems to lead to a segmentation fault
# degree=3


[sample]
mateFile = {snakemake.input.tumor_bam}
inputFormat = BAM
mateOrientation = 0

# If you use sorted SAM or BAM, please set "mateOrientation=0";
# then FREEC will not try to detect pairs with normal orientation and insert size.
# Instead, it will keep all pairs from the input file


[control]
mateFile = {snakemake.input.normal_bam}
inputFormat = BAM
mateOrientation = 0

# If you use sorted SAM or BAM, please set "mateOrientation=0";
# then FREEC will not try to detect pairs with normal orientation and insert size.
# Instead, it will keep all pairs from the input file

[BAF]
# no parameters to be set

[target]
# no parameters to be set
EOF

# Run Control-Freec -------------------------------------------------------------------------------

freec --conf $TMPDIR/tmp.d/freec.ini {extra_arguments}

# Control-FREEC names its files after the tumor BAM; only the ratios are declared
mv $output_dir/$(basename {snakemake.input.tumor_bam})_ratio.txt {snakemake.output.ratio}
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Control-FreeC: Snakemake wrapper.py
"""

import fnmatch

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

cf_config = snakemake.config["step_config"]["somatic_wgs_cnv_calling"]["control_freec"]

# Set parameters which depend on each other
defaults = {
    "contamination": "contaminationAdjustment = FALSE",
    "mappability": "",
    "window_size": ""
}

if "contamination" in cf_config and cf_config["contamination"] != "":
    defaults["contamination"] = "contaminationAdjustment = TRUE\ncontamination = {}".format(cf_config["contamination"])

if "path_mappability" in cf_config and cf_config["path_mappability"] != "":
    defaults["mappability"] = "gemMappabilityFile = {}".format(cf_config["path_mappability"])
    if "uniqueMatch" in cf_config and cf_config["uniqueMatch"]:
        defaults["mappability"] += "\n" + "uniqueMatch = TRUE"

if "coefficientOfVariation" in cf_config:
    defaults["window_size"] = "coefficientOfVariation = {}".format(cf_config["coefficientOfVariation"])
if "window_size" in cf_config and cf_config["window_size"] > 0:
    defaults["window_size"] = "window = {}".format(cf_config["window_size"])
    if "step" in cf_config and cf_config["step"] > 0:
        defaults["window_size"] += "\n" + "step = {}".format(cf_config["step"])

# Chromosome lengths without ignored chromosomes
lines = []
with open(snakemake.config["static_data_config"]["reference"]["path"] + ".fai", "rt") as f:
    for line in f:
        found = False
        for ignored in cf_config["ignore_chroms"]:
            if fnmatch.fnmatch(line[:line.index("\t")], ignored):
                found = True
                break
        if not found:
            lines.append(line.rstrip())
lines = "\n".join(lines)

print("DEBUG- snakemake.output = {}".format(snakemake.output), file=sys.stderr)

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Write out information about conda installation --------------------------------------------------

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Setup auto-cleaned TMPDIR -----------------------------------------------------------------------

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# write out configuration -------------------------------------------------------------------------
export output_prefix=$(basename {snakemake.output.ratio} .ratio.txt)
export output_dir=$(dirname {snakemake.output.ratio})

mkdir -p $TMPDIR/tmp.d

cat << __EOF > $TMPDIR/tmp.d/chrLen.fa.fai
{lines}
__EOF

cat <<EOF >$TMPDIR/tmp.d/freec.ini
# http://boevalab.com/FREEC/tutorial.html
[general]

# Part 1: user control of algorithm

breakPointThreshold = {cf_config[breakPointThreshold]}

chrLenFile = $TMPDIR/tmp.d/chrLen.fa.fai

{defaults[contamination]}

{defaults[mappability]}
minMappabilityPerWindow = {cf_config[minMappabilityPerWindow]}

minCNAlength = {cf_config[minCNAlength]}

{defaults[window_size]}

minExpectedGC = {cf_config[minExpectedGC]}
maxExpectedGC = {cf_config[maxExpectedGC]}

minimalSubclonePresence = {cf_config[minimalSubclonePresence]}

readCountThreshold = {cf_config[readCountThreshold]}

telocentromeric = {cf_config[telocentromeric]}


# Part 2: Hard-coded defaults

ploidy = 2

# this makes CNV predictions on Gonosomes unreliable
sex=XX

breakPointType = 2

forceGCcontentNormalization = 0

intercept = 0

noisyData = FALSE

printNA = TRUE

# this sets the degree of polynomial; currently the step of estimating it seems to lead to a segmentation fault
degree = 3

# chrFiles & GCcontentProfile not defined, GC correction using control sample
# chrFiles = <missing>
# GCcontentProfile = <missing>

# Part 3: Resources & output
maxThreads = 16
outputDir = $output_dir

BedGraphOutput = TRUE

bedtools = bedtools
sambamba = sambamba
samtools = samtools

SambambaThreards = 16

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

cp $TMPDIR/tmp.d/freec.ini $output_dir
freec --conf $TMPDIR/tmp.d/freec.ini

pushd $output_dir

for suffix in bam_control.cpn bam_CNVs bam_info.txt bam_ratio.BedGraph bam_ratio.txt bam_sample.cpn bam_subclones.txt; do
    new_suffix=$(echo ${{suffix}} | sed 's/bam_//')
    mv *.${{suffix}} $output_prefix.${{new_suffix}}
done

for f in *; do
    test -f $f && md5sum $f >$f.md5
done
popd

"""
)

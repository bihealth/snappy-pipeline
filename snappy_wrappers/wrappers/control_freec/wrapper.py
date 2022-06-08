# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Control-FreeC: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

window_str = ""
w = snakemake.config["step_config"]["somatic_wgs_cnv_calling"]["control_freec"]["window_size"]
if w >= 0:
    window_str = "window = {}".format(w)

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

cat <<EOF >$TMPDIR/tmp.d/freec.ini
# http://boevalab.com/FREEC/tutorial.html
[general]

maxThreads = 16
outputDir = $output_dir

## path to sambamba (faster BAM file reading)
sambamba = sambamba

chrLenFile = {snakemake.config[step_config][somatic_wgs_cnv_calling][control_freec][path_chrlenfile]}
ploidy = 2

breakPointThreshold = .8
coefficientOfVariation = 0.05
{window_str}

# information from pathology review is that samples are at least 60% tumor tissue
contaminationAdjustment = TRUE
contamination = 0.4
minimalSubclonePresence = 0.2

numberOfProcesses = 4

$(if [[ "{snakemake.config[step_config][somatic_wgs_cnv_calling][control_freec][path_mappability_enabled]}" == True ]]; then
  echo gemMappabilityFile = {snakemake.config[step_config][somatic_wgs_cnv_calling][control_freec][path_mappability]};
fi)

uniqueMatch = TRUE
forceGCcontentNormalization = 0

minCNAlength = 2

BedGraphOutput=TRUE

# this makes CNV predictions on Gonosomes unreliable
sex=XX

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

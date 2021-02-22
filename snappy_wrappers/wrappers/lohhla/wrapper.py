# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for LOHHLA: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt"

shell.executable("/bin/bash")


shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/out

mkdir -p $TMPDIR/bams
# need to link bams into the same directory
ln -sr {snakemake.input.normal_bam} -t $TMPDIR/bams
ln -sr {snakemake.input.normal_bai} -t $TMPDIR/bams
ln -sr {snakemake.input.tumor_bam} -t $TMPDIR/bams
ln -sr {snakemake.input.tumor_bai} -t $TMPDIR/bams

for i in $TMPDIR/bams/*; do
    mv $i $(echo $i | sed 's/bwa.//')
done

normal=$(echo $TMPDIR/bams/*-N1-*.bam)

lohhla_script=$(which lohhla)

#Rscript LOHHLAscript.R \
#Rscript lohhla \

jellyfish -h

Rscript $lohhla_script \
    --BAMDir $TMPDIR/bams/  \
    --HLAexonLoc /fast/groups/cubi/projects/biotools/OptiType/data/hla.dat \
    --HLAfastaLoc /fast/groups/cubi/projects/biotools/OptiType/data/hla_reference_dna.fasta \
    --cleanUp FALSE \
    --fishingStep FALSE \
    --gatkDir /fast/groups/cubi/projects/biotools/picard-tools-1.119 \
    --hlaPath $(realpath {snakemake.input.hla} ) \
    --mappingStep TRUE \
    --minCoverageFilter 10 \
    --normalBAMfile $normal \
    --novoDir $(dirname $(realpath $lohhla_script)) \
    --outputDir $(dirname $(realpath {snakemake.output.done} )) \
    --patientId example
    ##--CopyNumLoc $(realpath example-file/solutions.txt) \
    ##--normalBAMfile $(realpath {snakemake.input.normal_bam} ) \

#touch {snakemake.output.done}
"""
)

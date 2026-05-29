# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for LOHHLA: Snakemake wrapper.py"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Clemens Messerschmidt"

ShellWrapper(snakemake).run(
    r"""
# Setup subdirectories in the automatically created TMPDIR
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

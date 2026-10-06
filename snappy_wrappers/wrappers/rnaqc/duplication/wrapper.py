# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FeatureCounts: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

shell.executable("/bin/bash")

shell(
    r"""
set -euo pipefail
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# ----------------------------------------------------------------------------
# Inititalisation: paired & strandedness decisions
# ----------------------------------------------------------------------------
# Find out single or paired ended
n_pair=$(samtools view -f 0x1 {snakemake.input.bam} | head -n 1000 | wc -l || true)
if [[ $n_pair -eq 0 ]]; then
    paired=0
else
    paired=1
fi

# Find out strand
strand={args[strand]}

if [ ${{strand}} -eq -1 ]
then
    strand=$(cat {snakemake.input.decision})
fi

# ----------------------------------------------------------------------------
# Duplication (from RSeQC package
# ----------------------------------------------------------------------------

mkdir ${{TMPDIR}}/duplication
read_duplication.py -i {snakemake.input.bam} -o ${{TMPDIR}}/duplication/out

mv ${{TMPDIR}}/duplication/out.seq.DupRate.xls {snakemake.output.dupl_seq}
md5sum {snakemake.output.dupl_seq} > {snakemake.output.dupl_seq_md5}
mv ${{TMPDIR}}/duplication/out.pos.DupRate.xls {snakemake.output.dupl_pos}
md5sum {snakemake.output.dupl_pos} > {snakemake.output.dupl_pos_md5}
"""
)

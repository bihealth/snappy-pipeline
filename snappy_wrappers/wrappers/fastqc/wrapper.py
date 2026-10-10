# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for FastQC: Snakemake wrapper.py"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# The reads are inputs. If another task wrote them, the input is that task's .done file and
# params list the files.
reads = args.get("input", {})
reads_left = snakemake.input.get("reads_left") or reads["reads_left"]
reads_right = snakemake.input.get("reads_right") or reads.get("reads_right", [])
more_reads = [*reads_left, *reads_right]

ShellWrapper(snakemake).run(
    r"""

outdir={snakemake.output.html}
mkdir -p $outdir

_JAVA_OPTIONS="-Xms256m -Xmx512m -XX:CompressedClassSpaceSize=512m" \
fastqc \
    --noextract \
    -o $outdir \
    -t {args[num_threads]} \
    $(echo {more_reads} | tr ' ' '\n' | grep 'fastq.gz$\|fastq$\|sam$\|bam$')

pushd $outdir
pwd
ls -lh

# Index of the reports, which the Snakemake report opens
{{
    echo "<html><body><h1>FastQC</h1><ul>"
    for report in *_fastqc.html; do
        echo "<li><a href=\"$report\">$report</a></li>"
    done
    echo "</ul></body></html>"
}} > index.html
"""
)

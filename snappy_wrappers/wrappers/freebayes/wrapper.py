# -*- coding: utf-8 -*-
"""Candidate calling with freebayes over the libraries of an entity, in one genome window.

Genotyping is left to varlociraptor, so freebayes runs with ``--pooled-continuous``. With a
``targets`` input (the targets of the libraries' kit), the window is cut to the targets. A window
without targets is cut to its first base, so that its chunk has the header and the samples of the
other chunks and the gather can concatenate them.
"""

import os
import tempfile
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake


def read_bed(path):
    """Return the (contig, start, end) intervals of a BED file."""
    intervals = []
    with open(path, "rt") as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and not fields[0].startswith(("#", "track", "browser")):
                intervals.append((fields[0], int(fields[1]), int(fields[2])))
    return intervals


def intersect(windows, targets):
    """Return the parts of the targets that lie in the windows, in the order of the targets."""
    by_contig = {}
    for contig, start, end in windows:
        by_contig.setdefault(contig, []).append((start, end))
    result = []
    for contig, start, end in targets:
        for window_start, window_end in by_contig.get(contig, ()):
            begin, stop = max(start, window_start), min(end, window_end)
            if begin < stop:
                result.append((contig, begin, stop))
    return result


args = snakemake.params.args
windows = read_bed(snakemake.input.region)
if targets_path := snakemake.input.get("targets"):
    regions = intersect(windows, read_bed(targets_path))
    if not regions:
        contig, start, _ = windows[0]
        regions = [(contig, start, start + 1)]
else:
    regions = windows

with tempfile.NamedTemporaryFile("wt", suffix=".bed", delete=False) as f:
    targets = f.name
    for contig, start, end in regions:
        print(contig, start, end, sep="\t", file=f)

bams = " ".join(snakemake.input.bam)

try:
    ShellWrapper(snakemake).run(
        r"""
freebayes \
    --fasta-reference {snakemake.input.reference} \
    --targets {targets} \
    --pooled-continuous \
    --min-alternate-count {args[min_alternate_count]} \
    --min-alternate-fraction {args[min_alternate_fraction]} \
    {extra_arguments} \
    {bams} \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type b \
        --output {snakemake.output.bcf}
bcftools index --force {snakemake.output.bcf}
"""
    )
finally:
    os.unlink(targets)

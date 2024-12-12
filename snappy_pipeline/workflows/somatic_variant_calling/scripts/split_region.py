import csv
from bisect import bisect
from contextlib import redirect_stderr
from itertools import accumulate

from snappy_wrappers.tools.genome_windows import matches_any


def split_regions(
    contig_lengths: dict[str, int], n_chunks: int, offset: int = 1
) -> list[list[[str, int, int]]]:
    """Splits a given collection of contigs based on their length into `n_chunks` equal sized chunks"""
    contigs, lengths = zip(*contig_lengths.items())
    cumulative_lengths = list(accumulate(lengths, initial=0))
    chunk_size = cumulative_lengths[-1] // n_chunks

    all_entries = []
    cumulative_start = cumulative_end = 0
    for _ in range(n_chunks):
        entries = []
        idx = bisect(cumulative_lengths, cumulative_start) - 1
        accumulated_length = 0
        while accumulated_length < chunk_size and idx < len(contigs):
            contig = contigs[idx]
            available_length = cumulative_lengths[idx + 1] - cumulative_end
            remaining = chunk_size - accumulated_length
            cumulative_end = cumulative_start + min(available_length, remaining)

            start = cumulative_start - cumulative_lengths[idx]
            end = cumulative_end - cumulative_lengths[idx]
            entries.append([contig, start + offset, end + offset])

            diff = cumulative_end - cumulative_start
            accumulated_length += diff
            cumulative_start += diff
            idx = bisect(cumulative_lengths, cumulative_start) - 1

        assert sum(map(lambda x: x[2] - x[1], entries)) == chunk_size

        all_entries.append(entries)

    # make sure the last entry has any bases lost due to rounding
    all_entries[-1][-1][2] = lengths[-1]
    return all_entries


with open(snakemake.log.log, "wt") as log:
    with redirect_stderr(log):
        with open(snakemake.params.fai, "rt") as fai:
            csv_reader = csv.reader(fai, delimiter="\t")
            lengths = {
                contig: int(length)
                for contig, length, *_ in csv_reader
                if not matches_any(contig, snakemake.params.ignore_chroms)
            }
        regions = split_regions(lengths, len(snakemake.output.regions))
        for region, path in zip(regions, snakemake.output.regions):
            with open(path, "wt") as f:
                for chrom, start, end in region:
                    print(f"{chrom}:{start}-{end}", file=f)

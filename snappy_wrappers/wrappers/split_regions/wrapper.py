import csv
import fnmatch
import sys
import heapq
from collections import namedtuple
from contextlib import redirect_stderr, redirect_stdout
from math import ceil


def matches_any(query, patterns):
    for pattern in patterns:
        if fnmatch.fnmatch(query, pattern):
            return True
    return False


Interval = namedtuple("Interval", ["contig", "start", "end"])
type Region = list[Interval]


def partition_lengths_with_splits(
    contigs: dict[str, int], num_partitions: int, slack: float = 0.1
) -> list[Region]:
    if not contigs:
        return [[] for _ in range(num_partitions)]

    total_length = sum(contigs.values())
    target_partition_size = ceil(total_length / num_partitions)

    partitions: list[Region] = [[] for _ in range(num_partitions)]

    min_heap = [(0, i) for i in range(num_partitions)]
    heapq.heapify(min_heap)

    sorted_contigs = sorted(contigs.items(), key=lambda x: x[1], reverse=True)

    def add_interval(
        partition_index: int, contig: str, start: int, length: int, contig_length: int
    ):
        end_position = min(start + length, contig_length)
        partitions[partition_index].append(Interval(contig, start, end_position))

    for contig, length in sorted_contigs:
        current_start = 0
        remaining_length = length

        while remaining_length > 0:
            current_size, partition_idx = heapq.heappop(min_heap)

            space_to_fill = target_partition_size - current_size

            if space_to_fill <= 0 or (remaining_length - space_to_fill) < (slack * target_partition_size):
                chunk_to_add = remaining_length
            else:
                chunk_to_add = space_to_fill

            chunk_to_add = min(remaining_length, chunk_to_add)
            add_interval(partition_idx, contig, current_start, chunk_to_add, length)

            current_start += chunk_to_add
            remaining_length -= chunk_to_add

            heapq.heappush(min_heap, (current_size + chunk_to_add, partition_idx))

    return partitions


if __name__ == "__main__":
    if snakemake := locals().get("snakemake", None):
        log = lambda: open(snakemake.log.log, "wt")  # noqa: E731
        fai_path = snakemake.input.fai
        if args := getattr(snakemake.params, "args", None):
            ignore_chroms = args["ignore_chroms"]
            padding = args["padding"]
        else:
            ignore_chroms = snakemake.params.ignore_chroms
            padding = snakemake.params.padding
        output_regions = snakemake.output.regions
    else:
        log = lambda: sys.stderr  # noqa: E731
        fai_path = sys.argv[1]
        ignore_chroms = sys.argv[2].split(",")
        padding = int(sys.argv[3])
        output_regions = sys.argv[4:]

    with log() as log, redirect_stderr(log), redirect_stdout(log):
        with open(fai_path, "rt") as fai:
            csv_reader = csv.reader(fai, delimiter="\t")
            contigs = {
                contig: int(length)
                for contig, length, *_ in csv_reader
                if not matches_any(contig, ignore_chroms)
            }
        num_regions = len(output_regions)
        regions = partition_lengths_with_splits(contigs, num_regions)

        for i, region in enumerate(regions):
            size = sum(map(lambda r: r.end - r.start, region))
            print(f"Region {i}, size {size}: {region}")

        # Verify that all intervals for each contig sum up to the contig length:
        for contig, length in contigs.items():
            relevant_intervals = [
                interval for region in regions for interval in region if interval.contig == contig
            ]
            total_length = sum(interval.end - interval.start for interval in relevant_intervals)
            assert total_length == length, (
                f"Total length of intervals for contig {contig} is {total_length}, expected {length}"
            )

        for region, path in zip(regions, output_regions):
            with open(path, "wt") as f:
                for contig, start, end in region:
                    contig_length = contigs[contig]
                    # BED format is 0 based, end exclusive
                    padded_start = max(start - padding, 0)
                    padded_end = min(end + padding + 1, contig_length)
                    print(f"{contig}\t{padded_start}\t{padded_end}", file=f)

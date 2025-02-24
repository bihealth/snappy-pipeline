import csv
import fnmatch
import sys
from collections import namedtuple
from contextlib import redirect_stderr, redirect_stdout

from math import ceil


def matches_any(query, patterns):
    for pattern in patterns:
        if fnmatch.fnmatch(query, pattern):
            return True
    return False


type Partition = list[tuple[str, int]]
Interval = namedtuple("Interval", ["contig", "start", "end"])
type Region = list[Interval]


def partition_lengths_with_splits(
    contigs: dict[str, int], num_partitions: int, slack: float = 0.1
) -> list[Region]:
    total_length = sum(contigs.values())
    target_partition_size = ceil(total_length / num_partitions)
    slack_threshold = slack * target_partition_size

    partitions: list[Region] = [[] for _ in range(num_partitions)]
    partition_sizes = [0] * num_partitions

    sorted_contigs = list(sorted(contigs.items(), key=lambda x: x[1], reverse=True))

    def calulate_end_position(start: int, length: int, contig_length: int) -> int:
        return min(start + length, contig_length)

    def add_interval(
        partition_index: int, contig: str, start: int, length: int, contig_length: int
    ):
        end_position = calulate_end_position(start, length, contig_length)
        partitions[partition_index].append(Interval(contig, start, end_position))
        partition_sizes[partition_index] += length

    for contig, length in sorted_contigs:
        smallest = partition_sizes.index(min(partition_sizes))
        remaining_space = target_partition_size - partition_sizes[smallest]

        # Start position within the current contig
        current_start = 0

        if partition_sizes[smallest] + length > target_partition_size:
            if remaining_space > slack_threshold:
                add_interval(smallest, contig, current_start, remaining_space, length)
                current_start += remaining_space
                remaining_length = length - remaining_space

                while remaining_length > 0:
                    smallest = partition_sizes.index(min(partition_sizes))
                    remaining_space = target_partition_size - partition_sizes[smallest]

                    if remaining_space > 0:
                        space_to_add = min(remaining_length, remaining_space)
                        add_interval(smallest, contig, current_start, space_to_add, length)
                        current_start += space_to_add
                        remaining_length -= space_to_add
                    else:
                        break
            else:
                next_group_index = (smallest + 1) % num_partitions
                add_interval(next_group_index, contig, current_start, length, length)
        else:
            add_interval(smallest, contig, current_start, length, length)

    return partitions


if __name__ == "__main__":
    if snakemake := locals()["snakemake"]:
        log = lambda: open(snakemake.log.log, "wt")  # noqa: E731
        fai_path = snakemake.input.fai
        ignore_chroms = snakemake.params.ignore_chroms
        padding = snakemake.params.padding
        output_regions = snakemake.output.regions
    else:
        log = lambda: sys.stderr  # noqa: E731
        fai_path = sys.argv[1]
        ignore_chroms = sys.argv[3].split(",")
        padding = int(sys.argv[4])
        output_regions = sys.argv[5:]

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
                    padded_start = max(start - padding, 0) + 1
                    padded_end = min(end + padding + 1, contig_length)
                    print(f"{contig}:{padded_start}-{padded_end}", file=f)

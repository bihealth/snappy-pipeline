import csv
import re
import sys
import heapq
from collections import namedtuple
from contextlib import redirect_stderr, redirect_stdout
from math import ceil


def matches_any(query, patterns):
    """Whether any regular expression matches ``query`` (``re.search``, i.e. unanchored)."""
    for pattern in patterns:
        if re.search(pattern, query):
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

            if space_to_fill <= 0 or (remaining_length - space_to_fill) < (
                slack * target_partition_size
            ):
                chunk_to_add = remaining_length
            else:
                chunk_to_add = space_to_fill

            chunk_to_add = min(remaining_length, chunk_to_add)
            add_interval(partition_idx, contig, current_start, chunk_to_add, length)

            current_start += chunk_to_add
            remaining_length -= chunk_to_add

            heapq.heappush(min_heap, (current_size + chunk_to_add, partition_idx))

    return partitions


def read_targets(paths, contigs: dict[str, int]) -> dict[str, list[tuple[int, int]]]:
    """Return the merged target intervals of the BED files per contig, clipped to the contigs."""
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for path in paths:
        with open(path, "rt") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3 or fields[0].startswith(("#", "track", "browser")):
                    continue
                contig, start, end = fields[0], int(fields[1]), int(fields[2])
                if contig in contigs:
                    by_contig.setdefault(contig, []).append(
                        (max(start, 0), min(end, contigs[contig]))
                    )
    merged = {}
    for contig, intervals in by_contig.items():
        intervals.sort()
        result = []
        for start, end in intervals:
            if result and start <= result[-1][1]:
                result[-1] = (result[-1][0], max(result[-1][1], end))
            elif start < end:
                result.append((start, end))
        merged[contig] = result
    return merged


def partition_by_targets(
    contigs: dict[str, int], targets: dict[str, list[tuple[int, int]]], num_partitions: int
) -> list[Region]:
    """Partition the contigs into windows of equal target bases.

    The windows follow the contigs in order and cover them whole, so a window is a stretch of
    the genome (the end of a contig and the start of the next, if a cut lies between them);
    cuts fall inside a target where its bases are shared.
    """
    total = sum(end - start for intervals in targets.values() for start, end in intervals)
    if total == 0:
        return partition_lengths_with_splits(contigs, num_partitions)
    quota = total / num_partitions
    partitions: list[Region] = [[] for _ in range(num_partitions)]
    index = 0  # the window being filled
    done = 0  # target bases in the windows before the current one and in the current one
    for contig, length in contigs.items():
        window_start = 0
        for start, end in targets.get(contig, []):
            position = start
            while position < end and index < num_partitions - 1:
                needed = quota * (index + 1) - done
                if end - position < needed:
                    break
                cut = position + ceil(needed)
                partitions[index].append(Interval(contig, window_start, cut))
                done += cut - position
                position = cut
                window_start = cut
                index += 1
            done += end - position
        if window_start < length:
            partitions[index].append(Interval(contig, window_start, length))
    return partitions


def cut_windows(contigs: dict[str, int], window_size: int) -> list[Interval]:
    """Cut each contig into windows of ``window_size`` bases from its start; its last window is shorter."""
    return [
        Interval(contig, start, min(start + window_size, length))
        for contig, length in contigs.items()
        for start in range(0, length, window_size)
    ]


def partition_windows(windows: list[Interval], num_partitions: int) -> list[Region]:
    """Distribute the windows in order over the partitions, about equal bases and at least one window each."""
    if len(windows) < num_partitions:
        raise ValueError(
            f"{len(windows)} windows for {num_partitions} parts: lower the number of jobs or the window size"
        )
    total = sum(window.end - window.start for window in windows)
    partitions: list[Region] = [[] for _ in range(num_partitions)]
    index = 0  # the partition being filled
    done = 0  # bases in the partitions up to the current one
    for i, window in enumerate(windows):
        # Move on once this partition has its share or the windows left are one per later partition
        if partitions[index] and (
            done * num_partitions >= total * (index + 1)
            or len(windows) - i == num_partitions - index - 1
        ):
            index += 1
        partitions[index].append(window)
        done += window.end - window.start
    return partitions


def main(fai_path, ignore_contigs, padding, output_regions, target_paths=(), window_size=None):
    with open(fai_path, "rt") as fai:
        csv_reader = csv.reader(fai, delimiter="\t")
        contigs = {
            contig: int(length)
            for contig, length, *_ in csv_reader
            if not matches_any(contig, ignore_contigs)
        }
    num_regions = len(output_regions)
    if window_size:
        regions = partition_windows(cut_windows(contigs, window_size), num_regions)
    elif target_paths:
        regions = partition_by_targets(contigs, read_targets(target_paths, contigs), num_regions)
    else:
        regions = partition_lengths_with_splits(contigs, num_regions)

    for i, region in enumerate(regions):
        size = sum(map(lambda r: r.end - r.start, region))
        print(f"Region {i}, size {size}: {region}")

    # Verify that all intervals for each contig sum up to the contig length:
    covered = {contig: 0 for contig in contigs.keys()}
    for region in regions:
        for interval in region:
            covered[interval.contig] += interval.end - interval.start
    for contig, length in contigs.items():
        assert covered[contig] == length, (
            f"Total length of intervals for contig {contig} is {covered[contig]}, expected {length}"
        )

    for region, path in zip(regions, output_regions):
        with open(path, "wt") as f:
            for contig, start, end in region:
                contig_length = contigs[contig]
                # BED format is 0 based, end exclusive
                padded_start = max(start - padding, 0)
                padded_end = min(end + padding, contig_length)
                print(f"{contig}\t{padded_start}\t{padded_end}", file=f)


if __name__ == "__main__":
    if snakemake := locals().get("snakemake", None):
        from snappy_wrappers.snappy_wrapper import PythonWrapper

        args = snakemake.params.args
        PythonWrapper(snakemake).run(
            lambda: main(
                snakemake.input.fai,
                args["ignore_contigs"],
                args["padding"],
                snakemake.output.regions,
                list(snakemake.input.get("targets", [])),
                args.get("window_size"),
            )
        )
    else:
        with redirect_stderr(sys.stderr), redirect_stdout(sys.stderr):
            main(sys.argv[1], sys.argv[2].split(","), int(sys.argv[3]), sys.argv[4:])

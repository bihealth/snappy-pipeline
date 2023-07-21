#!/usr/bin/env python3
"""Generate json file for varfish import from DRAGEN json.
"""
from argparse import ArgumentParser
from collections import defaultdict
import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
import re
from typing import List, Optional, Tuple

COVERAGE_PATTERNS = {
    "idxstats": "**/*idxstats.txt",
    "mapping_metrics": "**/*mapping_metrics.csv",
    "fastqc_metrics": "**/*fastqc_metrics.csv",
    "overall_mean_cov": "**/*_overall_mean_cov.csv",
    "metrics": "**/*_coverage_metrics.csv",
    "fine_hist": "**/*_fine_hist.csv",
    "hist": "**/*_hist.csv",
    "cov_report": "**/*_cov_report.bed",
}

BLACKLIST = ["read_cov_report"]


@dataclass(frozen=True)
class Histogram:
    depth: int
    count: int

    @classmethod
    def from_text(cls, depth: str, count: str) -> "Histogram":
        return cls(depth=int(depth.strip("+")), count=int(count))


def get_region_id(name: str) -> Optional[str]:
    """Extract the region id from dragen qc filenames."""
    if "wgs" in name:
        return "wgs"
    if m := re.search(r"qc-coverage-region-(\d+)", name):
        return m.group(1)


def find_coverage_files(path: Path) -> dict:
    """Return coverage files ordered by region."""
    files = defaultdict(dict)
    general_files = {}
    for pattern_type, pattern in COVERAGE_PATTERNS.items():
        for matched_path in path.glob(pattern):
            if not any(b in matched_path.name for b in BLACKLIST):
                if region_id := get_region_id(matched_path.name):
                    files[region_id][pattern_type] = matched_path
                else:
                    general_files[pattern_type] = matched_path

    # add general files to all regions
    for region_id, region_files in files.items():
        files[region_id] = {**region_files, **general_files}
    return files


def get_mean_coverage(path: Path) -> float:
    """Load mean coverage."""

    with path.open() as f:
        _, mean_coverage = f.read().strip().split(", ")
    return float(mean_coverage)


def load_coverage_histogram(path: Path) -> list:
    with path.open() as f:
        reader = csv.DictReader(f)
        hist = [Histogram.from_text(row["Depth"], row["Overall"]) for row in reader]
    return hist


def load_bed(path: Path) -> list:
    with path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        data = list(reader)
    return data


def get_target_metrics(target_report_path: Path) -> dict:
    """Load report bed file and extract the target count and total target size."""
    report_data = load_bed(target_report_path)

    total_target_size = 0
    for entry in report_data:
        size = int(entry["end"]) - int(entry["start"])
        total_target_size += size

    return {
        "target_count": len(report_data),
        "total_target_size": total_target_size,
    }


def get_cumulative(
    fine_hist: List[Histogram], start: int = 0, end: int = 200, step: int = 10
) -> Tuple[dict, dict]:
    """Get total and relative cumulative histograms."""
    counts = {i: 0 for i in range(start, end, step)}

    total_count = 0
    for hist in fine_hist:
        total_count += hist.count
        for i in counts:
            if hist.depth >= i:
                counts[i] += hist.count

    ratios = {i: c / total_count for i, c in counts.items()}

    return counts, ratios


def compute_average_quality(fastqc_metrics: dict) -> float:
    total_qual = 0.0
    total_count = 0.0
    for (_, name), value in fastqc_metrics["READ MEAN QUALITY"].items():
        if m := re.search(r"Q(\d+) Reads", name):
            quality = float(m.group(1))
            value = float(value)
            total_qual += quality * value
            total_count += value
    return total_qual / total_count


def get_bamstats(metrics_path: Path, fastqc_path: Path) -> dict:
    mapping_metrics = defaultdict(dict)
    with metrics_path.open() as f:
        for line in f:
            line = line.strip()
            match line.count(","):
                case 4:
                    cat, region, name, value, ratio = line.split(",", 4)
                case 3:
                    cat, region, name, value = line.split(",", 3)
                    ratio = 0
                case _:
                    raise RuntimeError

            mapping_metrics[cat][(region, name)] = (value, ratio)

    fastqc_metrics = defaultdict(dict)
    with fastqc_path.open() as f:
        for line in f:
            cat, readname, name, value = line.split(",", 3)
            fastqc_metrics[cat][(readname, name)] = value

    summary = mapping_metrics["MAPPING/ALIGNING SUMMARY"]

    raw_total_sequences = int(summary[("", "Total input reads")][0])
    filtered_sequences = int(summary[("", "QC-failed reads")][0])
    sequences = raw_total_sequences - filtered_sequences

    mate_reads = int(summary[("", "Reads with mate sequenced")][0])
    reads_mapped = int(summary[("", "Mapped reads")][0])
    reads_mapped_and_paired = int(summary[("", "Paired reads (itself & mate mapped)")][0])
    reads_unmapped = int(summary[("", "Unmapped reads")][0])
    reads_properly_paired = int(summary[("", "Properly paired reads")][0])
    reads_duplicated = int(summary[("", "Number of duplicate marked reads")][0])
    reads_mq0 = int(summary[("", "Reads with MAPQ [ 0:10)")][0])  # TODO: this does not match 1:1
    reads_qc_failed = int(summary[("", "QC-failed reads")][0])

    non_primary_alignments = int(summary[("", "Secondary alignments")][0])

    read_length = int(float(summary[("", "Estimated read length")][0]))
    maximum_length = math.ceil(float(summary[("", "Estimated read length")][0]))

    average_quality = compute_average_quality(fastqc_metrics)

    bases_mapped = int(summary[("", "Mapped bases")][0])
    bases_softclipped = int(summary[("", "Soft-clipped bases")][0])
    bases_trimmed = int(summary[("", "Hard-clipped bases")][0])
    bases_q30 = int(summary[("", "Q30 bases")][0])
    bases_q30_wo_dups = int(summary[("", "Q30 bases (excl. dups & clipped bases)")][0])
    bases_q30_dups = bases_q30 - bases_q30_wo_dups
    bases_mapped = int(summary[("", "Mapped bases")][0])
    bases_mismatched_r1 = int(summary[("", "Mismatched bases R1")][0])
    bases_mismatched_r2 = int(summary[("", "Mismatched bases R2")][0])
    bases_mismatched = bases_mismatched_r1 + bases_mismatched_r2
    error_rate = bases_mismatched / bases_mapped

    mean_insert_size = float(summary[("", "Insert length: mean")][0])
    std_insert_size = float(summary[("", "Insert length: standard deviation")][0])

    pair_not_proper = int(summary[("", "Not properly paired reads (discordant)")][0])
    pair_diff_chrom = int(summary[("", "Paired reads mapped to different chromosomes")][0])

    return {
        "raw total sequences": raw_total_sequences,
        "filtered sequences": filtered_sequences,
        "sequences": sequences,
        "is sorted": 1,  # DRAGEN Sorting is enabled by default, should in the future check dragen replay
        "1st fragments": mate_reads / 2,
        "last fragments": mate_reads / 2,
        "reads mapped": reads_mapped,
        "reads mapped and paired": reads_mapped_and_paired,
        "reads_unmapped": reads_unmapped,
        "reads properly paired": reads_properly_paired,
        "reads pairsdfsed": mate_reads,
        "reads duplicated": reads_duplicated,
        "reads MQ0": reads_mq0,
        "reads QC failed": reads_qc_failed,
        "non-primary alignments": non_primary_alignments,
        "total length": read_length * mate_reads,
        "total first fragment length": read_length * mate_reads / 2,
        "total last fragment length": read_length * mate_reads / 2,
        "bases mapped": bases_mapped,
        "bases mapped (cigar)": bases_mapped - bases_softclipped,
        "bases trimmed": bases_trimmed,
        "bases duplicated": bases_q30_dups,
        "mismatches": bases_mismatched,
        "error rate": error_rate,
        "average length": read_length,
        "average first fragment length": read_length,
        "average last fragment length": read_length,
        "maximum length": maximum_length,
        "maximum first fragment length": maximum_length,
        "maximum last fragment length": maximum_length,
        "average quality": average_quality,
        "insert size average": mean_insert_size,
        "insert size standard deviation": std_insert_size,
        "inward oriented pairs": 0,  # does not exist on DRAGEN
        "outward oriented pairs": 0,  # does not exist on DRAGEN
        "pairs with other orientation": pair_not_proper,
        "pairs on different chromosomes": pair_diff_chrom,
        "percentage of properly paired reads (%)": float(summary[("", "Properly paired reads")][1]),
    }


def get_idxstats(idxstats_path: Path) -> dict:
    idxstats = {}
    with idxstats_path.open() as f:
        for line in f:
            contig, contig_length, mapped, unmapped = line.strip().split("\t", 3)
            idxstats[contig] = {
                "mapped": int(mapped),
                "unmapped": int(unmapped),
            }
    return idxstats


def create_varfish_json(sample_id: str, region_coverage_paths: dict) -> dict:
    """Create a json with the following keys, as needed for varfish import.

    sample_id:
      summary:
        mean_coverage: xx.xx
        target count: xxxxxx
        total target size: xxxxxxxx
      min_cov_base:
        "x": x.xx
        0-200 in 10 steps
      min_cov_target:
        "x": x.xx
        0-200 in 10 steps
    """
    mean_coverage = get_mean_coverage(region_coverage_paths["overall_mean_cov"])
    target_metrics = get_target_metrics(region_coverage_paths["cov_report"])
    _, relative = get_cumulative(load_coverage_histogram(region_coverage_paths["fine_hist"]))

    bamstats_fmt = get_bamstats(
        region_coverage_paths["mapping_metrics"], region_coverage_paths["fastqc_metrics"]
    )
    idxstats_fmt = get_idxstats(region_coverage_paths["idxstats"])

    relative_fmt = {str(k): v for k, v in relative.items()}

    output_data = {
        "summary": {
            "mean coverage": mean_coverage,
            **target_metrics,
        },
        "min_cov_base": relative_fmt,
        "min_cov_target": relative_fmt,
        "bamstats": bamstats_fmt,
        "idxstats": idxstats_fmt,
    }
    return {sample_id: output_data}


if __name__ == "__main__":
    parser = ArgumentParser(
        description="Extract DRAGEN QC Information into a format importable into Varfish."
    )
    parser.add_argument("--sample", required=True, help="sample id used for varfish import")
    parser.add_argument(
        "--region", default="3", help="region id pointing to relevant region qc for import"
    )
    parser.add_argument(
        "input_dir", type=Path, help="Input directory containing the file tree for a single case."
    )

    args = parser.parse_args()
    region_coverage_files = find_coverage_files(args.input_dir)[args.region]
    result_data = create_varfish_json(args.sample, region_coverage_files)

    print(json.dumps(result_data))

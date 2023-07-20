#!/usr/bin/env python3
"""Generate json file for varfish import from DRAGEN json.
"""
from argparse import ArgumentParser
from collections import defaultdict
import csv
from dataclasses import dataclass
import json
from pathlib import Path
import re
from typing import List, Optional, Tuple

COVERAGE_PATTERNS = {
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
    for pattern_type, pattern in COVERAGE_PATTERNS.items():
        for matched_path in path.glob(pattern):
            if not any(b in matched_path.name for b in BLACKLIST):
                if region_id := get_region_id(matched_path.name):
                    files[region_id][pattern_type] = matched_path
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

    relative_fmt = {str(k): v for k, v in relative.items()}

    output_data = {
        "summary": {
            "mean coverage": mean_coverage,
            **target_metrics,
        },
        "min_cov_base": relative_fmt,
        "min_cov_target": relative_fmt,
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

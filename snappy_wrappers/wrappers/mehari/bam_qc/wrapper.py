# -*- coding: utf-8 -*-
"""The VarFish table of the BAM quality control of a pedigree.

One row with a JSON object per library: the ``SN`` section of ``samtools stats``, ``samtools
idxstats``, and the target coverage of alfred (mean coverage, target size, the fraction of the
targets at the coverage thresholds). The quotes of the JSON are tripled, as varfish-cli reads it.
The param ``libraries`` names the library of each ``bamstats`` and ``idxstats`` input, in their order.
"""

import gzip
import json

#: The coverage thresholds of the target coverage fractions
COV_THRESHOLDS = [5, 15, 20, 30, 50, 100, 200, 300, 400, 500, 1000]


def number(text):
    """An ``int`` or ``float`` for a numeric ``samtools stats`` value, else the text."""
    for parse in (int, float):
        try:
            return parse(text)
        except ValueError:
            continue
    return text


def parse_bamstats(path):
    """The ``SN`` lines of ``samtools stats``: name -> value."""
    result = {}
    with open(path, "rt") as f:
        for line in f:
            if not line.startswith("SN\t"):
                continue
            _, name, value = line.rstrip("\n").split("\t")[:3]
            result[name.rstrip(":")] = number(value)
    return result


def parse_idxstats(path):
    """``samtools idxstats``: contig -> mapped and unmapped reads."""
    result = {}
    with open(path, "rt") as f:
        for line in f:
            contig, _, mapped, unmapped = line.rstrip("\n").split("\t")[:4]
            result[contig] = {"mapped": int(mapped), "unmapped": int(unmapped)}
    return result


def parse_alfred_qc(path):
    """The target coverage of alfred's QC JSON: sample -> summary and coverage fractions."""
    with gzip.open(path, "rt") as f:
        data = json.load(f)
    result = {}
    for sample in data["samples"]:
        columns = sample["summary"]["data"]["columns"]
        first_row = dict(zip(columns, sample["summary"]["data"]["rows"][0]))
        aligned_in_target = first_row["#AlignedBasesInBed"]
        target_size = first_row["#TotalBedBp"]
        for metric in sample["readGroups"][0]["metrics"]:
            if metric["id"] == "targetCoverage":
                coverages = metric["x"]["data"][0]["values"]
                fractions = metric["y"]["data"][0]["values"]
                break
        else:
            raise RuntimeError("Could not find targetCoverage metric in {}".format(path))
        min_cov_target = {
            cov: frac * 100 for cov, frac in zip(coverages, fractions) if cov in COV_THRESHOLDS
        }
        result[sample["id"]] = {
            "summary": {
                "mean coverage": aligned_in_target / target_size,
                "total target size": target_size,
            },
            "min_cov_target": min_cov_target,
        }
    return result


def bam_qc_table(libraries, bamstats, idxstats, alfred_qc):
    """The table: a header and one row whose ``bam_stats`` is the JSON of all libraries."""
    stats = {}
    for library, path in zip(libraries, bamstats, strict=True):
        stats.setdefault(library, {})["bamstats"] = parse_bamstats(path)
    for library, path in zip(libraries, idxstats, strict=True):
        stats.setdefault(library, {})["idxstats"] = parse_idxstats(path)
    for path in alfred_qc:
        for sample, coverage in parse_alfred_qc(path).items():
            stats.setdefault(sample, {}).update(coverage)
    bam_stats = json.dumps(stats, separators=(",", ":")).replace('"', '"""')
    return "case_id\tset_id\tbam_stats\n.\t.\t{}\n".format(bam_stats)


def main(libraries, bamstats, idxstats, alfred_qc, output):
    with gzip.open(output, "wt") as f:
        f.write(bam_qc_table(libraries, bamstats, idxstats, alfred_qc))


if __name__ == "__main__":
    if snakemake := locals().get("snakemake", None):
        from snappy_wrappers.snappy_wrapper import PythonWrapper

        PythonWrapper(snakemake).run(
            lambda: main(
                snakemake.params.args["libraries"],
                list(snakemake.input.bamstats),
                list(snakemake.input.idxstats),
                list(snakemake.input.get("alfred_qc", [])),
                snakemake.output.bam_qc,
            )
        )

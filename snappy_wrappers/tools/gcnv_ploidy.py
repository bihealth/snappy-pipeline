#!/usr/bin/env python3
"""Set the X and Y ploidy calls of GATK's DetermineGermlineContigPloidy to the sex of the PED files.

The calls are unreliable, which breaks JointGermlineCNVSegmentation further downstream
(https://github.com/broadinstitute/gatk/issues/8164).
Males get X 1 and Y 1, females X 2 and Y 0; ``--autosomes-diploid`` sets the autosomes to 2.
A changed call keeps the called ploidy in its quality, ``42.000<ploidy>``, as a marker.
The gCNV wrappers run this script with the Python of their environment.
"""

# The gCNV environment has Python 3.10, so this stays valid Python 3.7.

import argparse
from pathlib import Path

#: The ploidy of X and Y by the sex code of a PED file
PLOIDY = {"1": {"x": 1, "y": 1}, "2": {"x": 2, "y": 0}}


def read_sexes(ped_paths) -> dict:
    """The sex code of each sample of the PED files."""
    sexes = {}
    for path in ped_paths:
        with open(path, "rt") as inputf:
            for line in inputf:
                arr = line.strip().split("\t")
                if arr[4] not in PLOIDY:
                    raise ValueError(f"{path}: sample {arr[1]} has sex {arr[4]!r}, not 1 or 2")
                sexes[arr[1]] = arr[4]
    return sexes


def sex_chromosome(contig: str):
    """``x`` or ``y`` for the contigs X, chrX, Y and chrY, in any case; else ``None``."""
    name = contig.lower()
    name = name[3:] if name.startswith("chr") else name
    return name if name in ("x", "y") else None


def fixed_calls(lines, ploidy: dict, autosomes_diploid: bool) -> list:
    """The lines of a ``contig_ploidy.tsv`` with the ploidy of X and Y from ``ploidy``."""
    result = []
    for line in lines:
        arr = line.split("\t")
        if arr[0].startswith("@") or arr[0].lower() == "contig":
            expected = None  # header
        elif sex_chromosome(arr[0]):
            expected = ploidy[sex_chromosome(arr[0])]
        else:
            expected = 2 if autosomes_diploid else None
        if expected is not None:
            called = int(arr[1])
            if called != expected:
                arr[1] = str(expected)
                arr[2] = f"42.000{called}"  # magic marker
        result.append("\t".join(arr).strip())
    return result


def fix_ploidy_calls(ploidy_dir, sexes: dict, autosomes_diploid: bool) -> None:
    """Rewrite the ``contig_ploidy.tsv`` of each sample below ``ploidy_dir``."""
    for sample_dir in sorted((Path(ploidy_dir) / "ploidy-calls").glob("SAMPLE_*")):
        sample_name = (sample_dir / "sample_name.txt").read_text().strip()
        if sample_name not in sexes:
            raise ValueError(f"{sample_dir}: sample {sample_name} is in no PED file")
        path_call = sample_dir / "contig_ploidy.tsv"
        with path_call.open("rt") as inputf:
            lines = fixed_calls(inputf.readlines(), PLOIDY[sexes[sample_name]], autosomes_diploid)
        path_call.write_text("".join(line + "\n" for line in lines))


def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--autosomes-diploid", action="store_true", help="set the autosomes to ploidy 2"
    )
    parser.add_argument("ploidy_dir", help="output directory of DetermineGermlineContigPloidy")
    parser.add_argument("ped", nargs="+", help="PED files with the sex of every sample")
    args = parser.parse_args(argv)
    fix_ploidy_calls(args.ploidy_dir, read_sexes(args.ped), args.autosomes_diploid)


if __name__ == "__main__":
    main()

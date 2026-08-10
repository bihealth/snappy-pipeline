import json
import logging
import os
import re
import sys

from typing import Any

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"


CLASS_I_LOCII = ("A", "B", "C", "E", "F", "G")
CLASS_II_LOCII = (
    "DRA",
    "DRB1",
    "DRB3",
    "DRB4",
    "DRB5",
    "DQA1",
    "DQB1",
    "DOA",
    "DOB",
    "DMA",
    "DMB",
    "DPA1",
    "DPB1",
)
ALLELE_PATTERN = re.compile(
    r"^(?P<prefix>HLA-)?"
    r"(?P<locus>A|B|C|E|F|G|(DP|DQ)(A1|B1)|(DO|DM)[AB]|DRA|DRB[1345])"
    r"\*(?P<major>[0-9]+):(?P<minor>[0-9]+[GNQSCA]?)"
    r"(:(?P<extra>.+))?$"
)


def read_hla_values(hla_typing_file: str, mhc_class: str | None = None) -> set[str]:
    hla_types = {}

    with open(hla_typing_file, "rt") as f:
        calls: dict[str, Any] = json.load(f)

        for locus, alleles in calls.items():
            if locus in CLASS_I_LOCII:
                cls = "class_i"
            elif locus in CLASS_II_LOCII:
                cls = "class_ii"
            else:
                logging.warning(f"Unsupported locus {locus} found in file {fn}")
                continue
            if mhc_class and mhc_class != cls:
                continue

            if locus not in hla_types:
                hla_types[locus] = set()

            prefix = "HLA-" if cls == "class_i" else ""

            for allele in alleles:
                m = ALLELE_PATTERN.match(allele)
                if not m:
                    continue
                assert locus == m.group("locus"), f"Allele {allele} not in locus {locus}"
                cleaned = prefix + locus + "*" + m.group("major") + ":" + m.group("minor")
                hla_types[locus].add(cleaned)

    results = set()
    for pair_a, pair_b in (("DPA1", "DPB1"), ("DQA1", "DQB1")):
        for a in hla_types.get(pair_a, set()):
            for b in hla_types.get(pair_b, set()):
                results.add(f"{a}-{b}")

    for locus in hla_types.keys():
        if locus in ("DPA1", "DPB1", "DQA1", "DQB1"):
            continue
        results |= hla_types[locus]

    return results

logging.basicConfig(
    format="%(asctime)s %(levelname)s: %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
)

if "snakemake" in globals():
    hla_files = snakemake.input
    hla_table = snakemake.output.hla_types
else:
    hla_files = sys.argv[2:]
    hla_table = sys.argv[1]

hla_types = {}
for fn in hla_files:
    results = read_hla_values(str(fn))
    for hla_type in results:
        if hla_type not in hla_types:
            hla_types[hla_type] = []
        hla_types[hla_type].append(os.path.basename(str(fn)))

with open(str(hla_table), "wt") as f:
    f.write("HLA Allele\tN\tFiles\n")
    for hla_type in sorted(list(hla_types.keys())):
        fns = hla_types[hla_type]
        f.write(hla_type + "\t" + str(len(fns)) + "\t" + ";".join(fns) + "\n")
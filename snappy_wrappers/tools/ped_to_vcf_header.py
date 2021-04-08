# -*- coding: utf-8 -*-
"""Add PED file content to VCF file header.

Usage::

    $ snappy-ped_to_vcf_header --ped-file PED --output TXT
"""

from __future__ import print_function

import argparse
from collections import OrderedDict, defaultdict, namedtuple
import os
import re
import sys

__author__ = "Oliver Stolpe <oliver.stolpe@bihealth.de>"


#: Translation scheme for PED attributes (sex, disease) to text
PED_TRANSLATE = OrderedDict(
    [
        ("Sex", OrderedDict([("0", "Unknown"), ("1", "Male"), ("2", "Female")])),
        ("Disease", OrderedDict([("0", "Unknown"), ("1", "Unaffected"), ("2", "Affected")])),
    ]
)

#: Template VCF header string for PED attributes (sex, disease)
TPL_META = "##META=<ID={id},Type=String,Number=1,Values={values}>"
#: Template VCF header string for samples in PED
TPL_SAMPLE = "##SAMPLE=<ID={id},Sex={sex},Disease={disease}>"
#: Template VCF header string for pedigree structure
TPL_PEDIGREE = "##PEDIGREE=<ID={id},Family={family},Father={father},Mother={mother}>"

#: Donor representation
Donor = namedtuple("Donor", ["family", "id", "father", "mother", "sex", "disease"])


def parse_ped(ped_file):
    """Parse a given PED file and yield each line as a Donor."""
    for line in ped_file.readlines():
        line = re.split("\s+", line.rstrip())[:6]

        if line[0].startswith("#"):
            continue

        if not len(line) == 6:
            raise Exception("PED file not complete.")

        yield Donor(*line)


def ped_vcf_header(donors):
    """Return VCF header string given donors."""
    snippet = []
    families = defaultdict(list)

    for key, value in PED_TRANSLATE.items():
        snippet.append(TPL_META.format(id=key, values="[{}]".format(", ".join(value.values()))))

    for donor in donors:
        families[donor.family].append(donor)
        snippet.append(
            TPL_SAMPLE.format(
                id=donor.id,
                sex=PED_TRANSLATE["Sex"][donor.sex],
                disease=PED_TRANSLATE["Disease"][donor.disease],
            )
        )

    for _, members in families.items():
        for member in members:
            if len(members) == 1 or not member.father == "0" or not member.mother == "0":
                snippet.append(
                    TPL_PEDIGREE.format(
                        id=member.id,
                        family=member.family,
                        father=member.father,
                        mother=member.mother,
                    )
                )

    return "\n".join(snippet)


def write_header_snippet(header, output):
    """Open a text file (create folders if necessary) and write content."""
    if os.path.dirname(output) and not os.path.exists(os.path.dirname(output)):
        os.makedirs(os.path.dirname(output))

    with open(output, "wt") as fh:
        fh.write("{}\n".format(header))


def run(args):
    """Program entry point after parsing the command line."""
    donors = parse_ped(args.ped_file)
    header = ped_vcf_header(donors)
    write_header_snippet(header, args.output)


def main(argv=None):
    """Program entry point for parsing the command line."""
    parser = argparse.ArgumentParser(
        description=("Parse PED file and transform pedigree information into" "VCF header format")
    )
    parser.add_argument(
        "--ped-file",
        type=argparse.FileType("rt"),
        required=True,
        help="PED file that contains the pedigree information",
    )
    parser.add_argument(
        "--output", type=str, required=True, help="File with PED information as VCF header snippet"
    )

    args = parser.parse_args(argv)
    run(args)


if __name__ == "__main__":
    sys.exit(main())

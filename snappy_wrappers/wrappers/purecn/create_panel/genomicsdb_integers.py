"""Write the numbers of a GenomicsDB workspace that GATK writes as strings as integers.

The genomicsdb R package of the PureCN image (0.0.4) cannot read the workspace otherwise.
"""

import argparse
import json
import os

#: The fields by file of the workspace, and the list in that file whose entries have them
FIELDS = {
    "vidmap.json": ("contigs", ("tiledb_column_offset", "length")),
    "callset.json": ("callsets", ("row_idx", "idx_in_file")),
}


def write_integers(workspace: str) -> None:
    for name, (entries, fields) in FIELDS.items():
        path = os.path.join(workspace, name)
        with open(path, "rt") as f:
            content = json.load(f)
        for entry in content[entries]:
            for field in fields:
                entry[field] = int(entry[field])
        with open(path, "wt") as f:
            json.dump(content, f, indent=2)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workspace", help="The GenomicsDB workspace directory")
    write_integers(parser.parse_args().workspace)

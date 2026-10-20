import gzip
import re

from snappy_wrappers.snappy_wrapper import PythonWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

#: pVACsplice's junction step needs Ensembl gene ids without their version
PATTERN = re.compile(r'^(.+)gene_id "(ENSG[0-9]{11})\.[0-9]+([^"]*)"(.*)')


def main():
    features = str(snakemake.input.features)
    opener = gzip.open if features.endswith(".gz") else open
    with opener(features, "rt") as f, open(snakemake.output.gtf, "wt") as g:
        for line in f:
            line = line.strip()
            if m := PATTERN.match(line):
                line = f'{m.group(1)}gene_id "{m.group(2)}{m.group(3)}"{m.group(4)}'
            g.write(line + "\n")


PythonWrapper(snakemake).run(main)

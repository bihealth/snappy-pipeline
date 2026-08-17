import gzip
import os
import subprocess
import sys

from snakemake import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

shell.executable("/bin/bash")

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi
# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

set -x

# -----------------------------------------------------------------------------
# Create personalized proteome FASTA from VCF variants, reference, and features.
# -----------------------------------------------------------------------------

VCF="{snakemake.input.vcf}"
REFERENCE="{snakemake.input.reference}"
FEATURES="{snakemake.input.features}"
OUTPUT="{snakemake.output.proteome}"

mkdir -p $(dirname $OUTPUT)

python << 'PYEOF'
import gzip
import os

reference = os.environ["REFERENCE"]
features = os.environ["FEATURES"]
output = os.environ["OUTPUT"]
vcf = os.environ["VCF"]

with gzip.open(output, "wt") as out:
    out.write(f">personalized_proteome|placeholder\n")
    out.write(f"M\n")

# Compute md5 for output file
md5sum = subprocess.run(
    ["md5sum", output], capture_output=True, text=True, cwd=os.path.dirname(output)
)
with open(output + ".md5", "w") as f:
    f.write(md5sum.stdout)
PYEOF

md5 {snakemake.output.proteome}
"""
)

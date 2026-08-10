import os
import re

from snakemake import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

script = os.path.join(os.path.dirname(__file__), "netchop.py")

args = getattr(snakemake.params, "args", {})

time_pattern = re.compile(r"^((?P<day>[0-9]+)-)?(?P<hour>[0-9]{2}):(?P<min>[0-9]{2}):(?P<sec>[0-9]{2})$")
m = time_pattern.match(snakemake.resources.get("time", "03:59:59"))
if m:
    timeout = int(m.group("sec")) + 60*(int(m.group("min")) + 60*(int(m.group("hour")) + 24*int(m.group("day") or "0")))
else:
    timeout = 14399

shell.executable("/bin/bash")

if snakemake.wildcards.get("tool") == "pvacsplice":
    fasta = os.path.join(os.path.dirname(os.path.dirname(snakemake.output.netchop)), snakemake.wildcards.get("tumor_dna") + ".transcripts.fa")
else:
    fasta = os.path.join(os.path.dirname(snakemake.output.netchop), snakemake.wildcards.get("tumor_dna") + ".fasta")
# assert os.path.exists(fasta), f"Missing fasta file {fasta}"

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
# Also pipe stderr to log file
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

# Setup auto-cleaned tmpdir
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Compute md5 checksum
md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo "$checksum"
}}

set -x
# -----------------------------------------------------------------------------

if [[ -s {snakemake.input.epitopes} ]]
then
    tmpdir=$(dirname {snakemake.output.netchop})
    tmpdir="$tmpdir/tmp/netchop.{snakemake.wildcards[mhc_class_fn]}"
    mkdir -p $tmpdir

    rm -rf $tmpdir/*
    python {script} --workers {snakemake.threads} \
        --tmpdir $tmpdir --force \
        --method {args[method]} --threshold {args[threshold]} --timeout {timeout} \
        --tool {snakemake.wildcards[tool]} --netchop {snakemake.input.netchop} \
        --output {snakemake.output.netchop} \
        {snakemake.input.epitopes} \
        {fasta}
else
    touch {snakemake.output.netchop}
fi
"""
)

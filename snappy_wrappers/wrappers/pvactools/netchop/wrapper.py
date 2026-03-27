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
    timeout = int(m.group("sec")) + 60*(int(m.group("min")) + 60*(int(m.group("hour")) + 24*int(m.group("day"))))
else:
    timeout = 14399

shell.executable("/bin/bash")

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

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

tmpdir=$(dirname {snakemake.output.epitopes})
tmpdir="$tmpdir/tmp"
rm -rf $tmpdir
mkdir $tmpdir

python {script} --workers {snakemake.threads} \
    --tmpdir $tmpdir \
    --method {args.method} --threshold {args.threshold} --timeout {timeout} \
    --netchop {snakemake.input.netchop} \
    --output {snakemake.output.epitopes} \
    {snakemake.input.vcf} \
    {snakemake.input.epitopes}
"""
)

"""CUBI+Snakemake wrapper code for scramble (cluster): Snakemake wrapper.py"""

from snakemake import shell

shell.executable("/bin/bash")


shell(
    r"""
set -x

# Pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Create out dir
mkdir -p $(dirname {snakemake.output.txt})

# Call tool
cluster_identifier {snakemake.input} > {snakemake.output.txt}
"""
)

# Compute MD5 sums of log
shell(
    r"""
md5sum {snakemake.log} > {snakemake.log}.md5
"""
)

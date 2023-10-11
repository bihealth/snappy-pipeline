"""CUBI+Snakemake wrapper code for scarHRD (non-conda package installation)
"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

lib_path = os.path.dirname(snakemake.input.done)

step = snakemake.config["pipeline_step"]["name"]
genome = snakemake.config["static_data_config"]["reference"]["path"]
length = snakemake.config["step_config"][step]["scarHRD"]["length"]
genome_name = snakemake.config["step_config"][step]["scarHRD"]["genome_name"]

chr_in_name = "TRUE" if snakemake.config["step_config"][step]["scarHRD"]["chr_prefix"] else "FALSE"
prefix = "chr" if snakemake.config["step_config"][step]["scarHRD"]["chr_prefix"] else ""
if genome_name == "grch37" or genome_name == "grch38":
    chromosomes = " ".join([prefix + str(x) for x in list(range(1, 23)) + ["X", "Y"]])
elif genome_name == "mouse":
    chromosomes = " ".join([prefix + str(x) for x in list(range(1, 21)) + ["X", "Y"]])
else:
    raise Exception("Invalid configuration")

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

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

export R_LIBS_USER="{lib_path}"
export VROOM_CONNECTION_SIZE=2000000000

export R_LIBS_PATH="{lib_path}"
export VROOM_CONNECTION_SIZE=2000000000

cat << __EOF | R --vanilla --slave
library("scarHRD")

tbl <- scar_score("{snakemake.input.seqz}", reference="{genome_name}", seqz=TRUE, chr.in.name={chr_in_name}, outputdir=dirname("{snakemake.output.scarHRD}"))
cat('{{\n', file="{snakemake.output.scarHRD}")
cat('    "HRD": ', tbl[1,1], ',\n', sep="", file="{snakemake.output.scarHRD}", append=TRUE)
cat('    "Telomeric AI": ', tbl[1,2], ',\n', sep="", file="{snakemake.output.scarHRD}", append=TRUE)
cat('    "LST": ', tbl[1,3], ',\n', sep="", file="{snakemake.output.scarHRD}", append=TRUE)
cat('    "HRD-sum": ', tbl[1,4], '\n', sep="", file="{snakemake.output.scarHRD}", append=TRUE)
cat('}}\n', file="{snakemake.output.scarHRD}", append=TRUE)

__EOF

pushd $(dirname {snakemake.output.scarHRD}) ; f=$(basename {snakemake.output.scarHRD}) ; md5sum $f > $f.md5 ; popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

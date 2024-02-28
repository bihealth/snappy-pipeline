"""CUBI+Snakemake wrapper code for sequenza (R part, post-processing)
"""

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snakemake import shell

from snappy_wrappers.tools.genome_windows import yield_contigs

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"


def config_to_r(x):
    if x is None:
        return "NULL"
    if isinstance(x, str):
        return f'"{x}"'
    if isinstance(x, bool):
        return "TRUE" if x else "FALSE"
    if isinstance(x, list):
        return "c({})".format(", ".join([config_to_r(xx) for xx in x]))
    if isinstance(x, dict):
        return "list({})".format(
            ", ".join(['"{}"={}'.format(k, config_to_r(v)) for k, v in x.items()])
        )
    return str(x)


step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["sequenza"]
genome = snakemake.config["static_data_config"]["reference"]["path"]

f = open(genome + ".fai", "rt")
contigs = config_to_r(list(yield_contigs(f, config.get("ignore_chroms"))))
f.close()

args_extract = config_to_r(dict(config["extra_args_extract"]))
args_fit = config_to_r(dict(config["extra_args_fit"]))

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

export R_LIBS_USER=$(dirname {snakemake.input.packages})
export VROOM_CONNECTION_SIZE=2000000000

R --vanilla --slave << __EOF
library(sequenza)

args <- list(file="{snakemake.input.seqz}", assembly="{config[assembly]}", chromosome.list={contigs})
args <- c(args, {args_extract})
seqz <- do.call(sequenza.extract, args=args)

args <- list(sequenza.extract=seqz, chromosome.list={contigs}, mc.cores=1)
args <- c(args, {args_fit})
CP <- do.call(sequenza.fit, args=args)

sequenza.results(sequenza.extract=seqz, cp.table=CP, sample.id="{snakemake.wildcards[library_name]}", out.dir=dirname("{snakemake.output.done}"))

__EOF

pushd $(dirname {snakemake.output.done}) ; fns=$(ls) ; for f in $fns ; do md5sum $f > $f.md5 ; done ; popd

touch {snakemake.output.done}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

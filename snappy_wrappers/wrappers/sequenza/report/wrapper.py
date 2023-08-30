"""CUBI+Snakemake wrapper code for sequenza (R part, post-processing)
"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

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

R --vanilla --slave << __EOF
library(sequenza)

seqz <- sequenza.extract("{snakemake.input.seqz}")
CP <- sequenza.fit(seqz)
sequenza.results(sequenza.extract=seqz, cp.table=CP, sample.id="{snakemake.params[sample_id]}", out.dir=dirname("{snakemake.output.done}"))

__EOF

pushd $(dirname {snakemake.output.done}) ; for f in $(ls) ; do md5sum $f > $f.md5 ; done ; popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

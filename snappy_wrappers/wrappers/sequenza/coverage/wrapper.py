"""CUBI+Snakemake wrapper code for sequenza (sequenza-utils, pileups)"""

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

args = getattr(snakemake.params, "args", {})

genome = args["reference"]
length = args["length"]

f = open(genome + ".fai", "rt")
contigs = " ".join(yield_contigs(f, args.get("ignore_chroms")))
f.close()

extra_arguments = " ".join(
    ["--{} {}".format(k, v) for k, v in args.get("extra_arguments", {}).items()]
)

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

sequenza-utils bam2seqz \
    -gc {snakemake.input.gc} --fasta {genome} \
    -n {snakemake.input.normal_bam} --tumor {snakemake.input.tumor_bam} \
    -C {contigs} {extra_arguments} \
    | sequenza-utils seqz_binning -s - \
        -w {length} -o {snakemake.output.seqz}

pushd $(dirname {snakemake.output.seqz}) ; f=$(basename {snakemake.output.seqz}) ; md5sum $f > $f.md5 ; popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

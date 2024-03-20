# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MANTIS: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Clemens Messerschmidt"

shell.executable("/bin/bash")


shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/out

# The following config is recommended by the authors for WES,
# but should also work for reasonable deep WGS according to them.
# https://github.com/OSU-SRLab/MANTIS/issues/25
#
# Only one thread is used, see https://github.com/OSU-SRLab/MANTIS/issues/57

mantis-msi2 \
    -t {snakemake.input.tumor_bam}  \
    -n {snakemake.input.normal_bam} \
    --genome {snakemake.config[static_data_config][reference][path]} \
    --bedfile {snakemake.config[step_config][somatic_msi_calling][loci_bed]} \
    --min-read-length 35 \
    --min-read-quality 20.0 \
    --min-locus-quality 25.0 \
    --min-locus-coverage 20 \
    --min-repeat-reads 1 \
    --threads 1 \
    -o $TMPDIR/out/$(basename {snakemake.output.result})

pushd $TMPDIR/out
for f in *; do
    md5sum $f > $f.md5
done
popd

mv $TMPDIR/out/* $(dirname {snakemake.output.result})
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

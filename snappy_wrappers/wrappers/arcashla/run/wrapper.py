# -*- coding: utf-8 -*-
"""Wrapper for actually running ASCAT"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

single = "--single" if args.get("single", False) else ""

extra_args = (
    f"--population {args['population']} --min_count {args['min_count']} --tolerance {args['tolerance']} "
    f"--max_iterations {args['max_iterations']} "
    f"--drop_threshold {args['drop_threshold']} "  # fails with --zygocity_threshold {args['zygocity_threshold']}" 
)
if "drop_iterations" in args:
    extra_args += f" --drop_iterations {args['drop_iterations']}"

if args.get("single", False):
    single = "--single"
    extra_args += f" --avg {args['avg']} --std {args['std']}"
else:
    single = ""

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

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

mkdir $TMPDIR/tmp
arcasHLA extract --verbose --threads {snakemake.threads} {single} \
    --temp $TMPDIR/tmp \
    --outdir $TMPDIR/extracted \
    {snakemake.input.bam}

rm -rf $TMPDIR/tmp

mkdir $TMPDIR/tmp
arcasHLA genotype --verbose --threads {snakemake.threads} \
    --outdir $TMPDIR/genotyped \
    --temp $TMPDIR/tmp \
    {extra_args} \
    $TMPDIR/extracted/*.fq.gz

cp $TMPDIR/genotyped/{args[mapper]}.genotype.json {snakemake.output.json}
pushd $(dirname {snakemake.output.json})
md5sum $(basename {snakemake.output.json}) >$(basename {snakemake.output.json}).md5
popd

# To make the pipeline happy.
# TODO: move all the step output to json.
md5sum {snakemake.log.log} > {snakemake.log.log}.md5
touch {snakemake.output.txt}
touch {snakemake.output.txt}.md5
"""
)

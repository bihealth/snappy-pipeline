# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py plot"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""

scatters = [
    getattr(snakemake.output, x)
    for x in filter(
        lambda x: x.startswith("scatter_chr") and not x.endswith("_md5"),
        [k for k, _ in snakemake.output._get_names()],
    )
]

ShellWrapper(snakemake).run(
    r"""
set -x

# -----------------------------------------------------------------------------

unset DISPLAY

if [[ -n "{snakemake.output.diagram}" ]]
then
    cnvkit.py diagram \
        --output {snakemake.output.diagram} \
        --segment {snakemake.input.cns} \
        {gender} {male} \
        --threshold {args[threshold]} --min-probes {args[min_probes]} \
        $(if [[ "{args[shift_xy]}" = "False" ]]; then \
            echo --no-shift-xy
        fi) \
        {snakemake.input.cnr}
else
    touch {snakemake.output.diagram}
fi

if [[ -n "{snakemake.output.scatter}" ]]
then
    cnvkit.py scatter \
        --output {snakemake.output.scatter} \
        --segment {snakemake.input.cns} \
        {gender} {male} \
        {snakemake.input.cnr}
else
    touch {snakemake.output.scatter}
fi

for scatter in {scatters}
do
    if [[ -n "$scatter" ]]
    then
        chrom=$(basename -s .png $scatter | sed -re "s/.*\.(chr([0-9]+|[XY]))$/\1/")
        cnvkit.py scatter \
            --output $scatter \
            --chromosome $chrom \
            --segment {snakemake.input.cns} \
            {gender} {male} \
            {snakemake.input.cnr}
    else
        touch $scatter
    fi
done
"""
)


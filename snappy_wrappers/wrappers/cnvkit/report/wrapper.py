# -*- coding: utf-8 -*-
"""Wrapper for cnvkit.py report"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})

gender = " --gender {}".format(args["gender"]) if args.get("gender", None) else ""
male = " --male-reference" if args.get("male_reference", False) else ""

drop_low_coverage = "--drop-low-coverage" if args.get("drop_low_coverage", False) else ""

input_target = getattr(snakemake.input, "target", "")
input_antitarget = getattr(snakemake.input, "antitarget", "")
input_cnr = getattr(snakemake.input, "cnr", "")
input_cns = getattr(snakemake.input, "cns", "")

# Each report tool is treated separately, as the arguments may not be in snakemake.params for all tools
if getattr(snakemake.output, "breaks", ""):
    if input_cnr != "" and input_cns != "":
        breaks = "cnvkit.py breaks --output {out} --min-probes {args[breaks][min_probes]}  {input_cnr} {input_cns}".format(
            input_cnr=input_cnr, input_cns=input_cns, args=args, out=snakemake.output.breaks
        )
    else:
        breaks = "touch {out}".format(out=snakemake.output.breaks)
else:
    breaks = ""

if getattr(snakemake.output, "genemetrics", ""):
    if input_cnr != "" and input_cns != "":
        genemetrics = """
cnvkit.py genemetrics \
    --output {out} \
    --segment {input_cns} \
    --min-probes {args[genemetrics][min_probes]} --threshold {args[genemetrics][threshold]} \
    {drop_low_coverage} {gender} {male} \
    --mean --median --mode --ttest --stdev --sem --mad --mse --iqr --bivar --ci --pi \
    --alpha {args[genemetrics][alpha]} --bootstrap {args[genemetrics][bootstrap]} \
    {input_cnr}
        """.format(
            input_cnr=input_cnr,
            input_cns=input_cns,
            args=args,
            out=snakemake.output.genemetrics,
            drop_low_coverage=drop_low_coverage,
            gender=gender,
            male=male,
        )
    else:
        genemetrics = "touch {out}".format(out=snakemake.output.genemetrics)
else:
    genemetrics = ""

if getattr(snakemake.output, "segmetrics", ""):
    if input_cnr != "" and input_cns != "":
        segmetrics = """
cnvkit.py segmetrics \
    --output {out} \
    --segment {input_cns} \
    --mean --median --mode --t-test --stdev --sem --mad --mse --iqr --bivar --ci --pi \
    --alpha {args[segmetrics][alpha]} --bootstrap {args[segmetrics][bootstrap]} \
    {smooth_bootstrap} {drop_low_coverage} \
    {input_cnr}
        """.format(
            input_cnr=input_cnr,
            input_cns=input_cns,
            args=args,
            out=snakemake.output.segmetrics,
            drop_low_coverage=drop_low_coverage,
            smooth_bootstrap="--smooth-bootstrap"
            if args["segmetrics"].get("smooth_bootstrap", False)
            else "",
        )
    else:
        segmetrics = "touch {out}".format(out=snakemake.output.segmetrics)
else:
    segmetrics = ""

if getattr(snakemake.output, "sex", ""):
    if input_cnr != "" or input_cns != "" or input_target != "" or input_antitarget != "":
        sex = """
        cnvkit.py sex \
            --output {out} {male} \
            {input_target} {input_antitarget} {input_cnr} {input_cns}
        """.format(
            input_cnr=input_cnr,
            input_cns=input_cns,
            input_target=input_target,
            input_antitarget=input_antitarget,
            args=args,
            out=snakemake.output.sex,
            male=male,
        )
    else:
        sex = "touch {out}".format(out=snakemake.output.sex)
else:
    sex = ""

if getattr(snakemake.output, "metrics", ""):
    if input_cnr != "" or input_cns != "" or input_target != "" or input_antitarget != "":
        metrics = """
cnvkit.py metrics \
    --output {out} \
    {drop_low_coverage} \
    {input_target} {input_antitarget} {input_cnr} {input_cns}
        """.format(
            input_cnr=input_cnr,
            input_target=input_target,
            input_antitarget=input_antitarget,
            input_cns=f"--segments {input_cns}" if input_cns else "",
            args=args,
            out=snakemake.output.metrics,
            drop_low_coverage=drop_low_coverage,
        )
    else:
        metrics = "touch {out}".format(out=snakemake.output.metrics)
else:
    metrics = ""

ShellWrapper(snakemake).run(
    r"""
set -x

{breaks}

{genemetrics}

{segmetrics}

{sex}

{metrics}
"""
)


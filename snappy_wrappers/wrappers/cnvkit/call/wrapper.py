# -*- coding: utf-8 -*-
"""Wrapper vor cnvkit.py call"""

import re

from snappy_wrappers.wrappers.cnvkit.cnvkit_wrapper import CnvkitWrapper

class CnvkitWrapperCall(CnvkitWrapper):
    PURITY_PATTERN = re.compile("^Purity: +([+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([EeDd][+-]?[0-9]+)?) *$")
    PLOIDY_PATTERN = re.compile("^Ploidy: +([+-]?([0-9]+(\.[0-9]*)?|\.[0-9]+)([EeDd][+-]?[0-9]+)?) *$")

    def preamble(self):
        if "purity" in self.snakemake.input:
            with open(self.snakemake.input.purity, "rt") as f:
                for line in f:
                    m = CnvkitWrapperCall.PURITY_PATTERN.match(line.strip())
                    if m:
                        self.purity = float(m.groups()[1])
                    else:
                        m = CnvkitWrapperCall.PLOIDY_PATTERN.match(line.strip())
                        if m:
                            self.ploidy = float(m.groups()[1])
        else:
            self.purity = self.snakemake.params.purity if "purity" in self.snakemake.params else None
            self.ploidy = self.snakemake.params.ploidy if "ploidy" in self.snakemake.params else None

        self.cmd = self.cmd.format(purity=self.purity, ploidy=self.ploidy)

if "variants" in snakemake.input:
    variants = r"""
        ---vcf {snakemake.input.variants} \
        {snakemake.params.sample_id} {snakemake.params.normal_id} \
        {snakemake.params.min_variant_depth} {snakemake.params.zygocity_freq}
    """.format(
        snakemake=snakemake,
    )
else:
    variants = ""

cmd = r"""
cnvkit.py call \
    -o {snakemake.output.calls} \
    --method {snakemake.params.method} --thresholds={snakemake.params.thresholds} \
    --filter {snakemake.params.filter} \
    {center} \
    {drop_low_coverage} \
    {sample_sex} {male_reference} \
    {variants} \
    {{purity}} {{ploidy}} \
    {snakemake.input.segments}
""".format(
    snakemake=snakemake,
    center=f"--center-at {snakemake.params.center_at}" if "center_at" in snakemake.params else f"--center {snakemake.params.center}",
    drop_low_coverage="--drop-low-coverage" if snakemake.params.drop_low_coverage else "",
    sample_sex=f"--sample-sex {snakemake.params.sample_sex}" if "sample_sex" in snakemake.params else "",
    male_reference="--male-reference" if snakemake.params.male_reference else "",
    variants=variants,
)

CnvkitWrapperCall(snakemake, cmd).run()

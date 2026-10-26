"""Compare the ngs-chew fingerprints of the task's libraries and plot the comparison"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

ShellWrapper(snakemake).run(
    r"""
# ngs-chew compare writes <prefix>.relatedness.txt, <prefix>.ibs0.txt and <prefix>.mask.txt
ngs-chew compare --output-prefix $TMPDIR/cohort {snakemake.input.npz}
mv $TMPDIR/cohort.relatedness.txt {snakemake.output.relatedness}
mv $TMPDIR/cohort.ibs0.txt {snakemake.output.ibs0}
mv $TMPDIR/cohort.mask.txt {snakemake.output.mask}
ngs-chew stats --output {snakemake.output.stats} {snakemake.input.npz}

ngs-chew plot-compare {snakemake.output.relatedness} {snakemake.output.relatedness_html}
ngs-chew plot-var-het {snakemake.output.stats} {snakemake.output.var_het_html}
"""
)

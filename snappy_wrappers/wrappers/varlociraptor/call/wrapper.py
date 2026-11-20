# -*- coding: utf-8 -*-
"""Call the variants of one chunk of the candidates of an entity against its scenario.

The scenario comes rendered as a param, so a changed scenario reruns the calls.
The wrapper writes it to the job's temporary directory and to the log.
With ``full_prior``, the prior weighs every event instead of only absent against present.
``varlociraptor genotype`` adds the genotypes that the consumers of germline and somatic calls expect.
The VCF keeps the observed samples only:
an unobserved sample of the scenario (an unsequenced relative, the normal of a tumor without one) has no data.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
observations = list(snakemake.input.observations)
assert len(observations) == len(args["obs_names"]), "one observation file per sample"
obs = " ".join("{}={}".format(name, path) for name, path in zip(args["obs_names"], observations))
samples = ",".join(args["obs_names"])
# The base class formats the command, so the braces of the scenario are doubled
scenario = args["scenario"].replace("{", "{{").replace("}", "}}")
full_prior = "--full-prior" if args["full_prior"] else ""

ShellWrapper(snakemake).run(
    r"""
cat > $TMPDIR/scenario.yaml << 'SCENARIO'
{scenario}
SCENARIO
cat $TMPDIR/scenario.yaml

varlociraptor call variants \
    {full_prior} \
    {extra_arguments} \
    generic \
    --obs {obs} \
    --scenario $TMPDIR/scenario.yaml \
    | varlociraptor genotype \
    | bcftools view \
        --samples {samples} \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type b \
        --output {snakemake.output.bcf}
bcftools index --force {snakemake.output.bcf}
"""
)

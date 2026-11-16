# -*- coding: utf-8 -*-
"""Call the variants of one genome window of an entity against its scenario.

The scenario comes rendered as a param, so a changed scenario reruns the calls; the wrapper
writes it to the job's temporary directory and to the log. ``varlociraptor genotype`` adds the
genotypes that the consumers of germline and somatic calls expect.
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
observations = list(snakemake.input.observations)
assert len(observations) == len(args["obs_names"]), "one observation file per sample"
obs = " ".join("{}={}".format(name, path) for name, path in zip(args["obs_names"], observations))
# The base class formats the command, so the braces of the scenario are doubled
scenario = args["scenario"].replace("{", "{{").replace("}", "}}")

ShellWrapper(snakemake).run(
    r"""
cat > $TMPDIR/scenario.yaml << 'SCENARIO'
{scenario}
SCENARIO
cat $TMPDIR/scenario.yaml

varlociraptor call variants \
    {extra_arguments} \
    generic \
    --obs {obs} \
    --scenario $TMPDIR/scenario.yaml \
    | varlociraptor genotype \
    | bcftools sort \
        --temp-dir $TMPDIR \
        --output-type b \
        --output {snakemake.output.bcf}
bcftools index --force {snakemake.output.bcf}
"""
)

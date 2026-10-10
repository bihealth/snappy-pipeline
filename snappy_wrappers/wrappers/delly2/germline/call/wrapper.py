from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
delly2_config = args["config"]

if exclude := snakemake.input.get("exclude", ""):
    exclude_str = f"--exclude {exclude}"
else:
    exclude_str = ""

ShellWrapper(snakemake).run(
    r"""
set -x

delly call \
    --map-qual {delly2_config[map_qual]} \
    --qual-tra {delly2_config[qual_tra]} \
    --geno-qual {delly2_config[geno_qual]} \
    --mad-cutoff {delly2_config[mad_cutoff]} \
    --genome {snakemake.input.reference} \
    --outfile {snakemake.output.bcf} \
    {exclude_str} \
    {snakemake.input.bam}

tabix -f {snakemake.output.bcf}
"""
)

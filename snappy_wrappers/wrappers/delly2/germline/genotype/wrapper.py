from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
delly2_config = args["config"]

if delly2_config["path_exclude_tsv"]:
    exclude_str = "--exclude %s" % delly2_config["path_exclude_tsv"]
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
    --vcffile {snakemake.input.bcf} \
    --genome {args[genome]} \
    --outfile {snakemake.output.bcf} \
    {exclude_str} \
    {snakemake.input.bam}

tabix -f {snakemake.output.bcf}
"""
)

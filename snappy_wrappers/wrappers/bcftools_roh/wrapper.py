from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

if path_targets := getattr(snakemake.input, "path_targets", ""):
    path_targets = f"--regions-file {path_targets}"
if path_af_file := getattr(snakemake.input, "path_af_file", ""):
    path_af_file = f"--AF-file {path_af_file}"

ignore_hormef = "--ignore-homref" if args.get("ignore_homref", False) else ""
skip_indels = "--skip-indels" if args.get("skip_indels", False) else ""
rec_rate = f"--rec_rate {args['rec_rate']}" if args.get("rec_rate", 0.0) > 0.0 else ""

ShellWrapper(snakemake).run(
    r"""
out={snakemake.output.txt}
raw_out=${{out%.regions.txt.gz}}.raw.txt.gz

bcftools roh \
    {path_targets} {path_af_file} \
    {ignore_homref} {skip_indels} {rec_rate} \
    --output $raw_out \
    --output-type srz \
    {snakemake.input.vcf}

# Cut out text and BED files.
(
    set +o pipefail
    zcat $raw_out \
    | head -n 3
    zcat $raw_out \
    | tail -n +4 \
    | egrep "^RG|^# RG"
) | bgzip -c \
> {snakemake.output.txt}
"""
)

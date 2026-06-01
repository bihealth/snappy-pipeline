# -*- coding: utf-8 -*-
"""Snakemake wrapper for hard-filtering the de novo results.

Apply all hard-filters with the following exceptions :

- don't remove variants with dbSNP IDs
- don't remove the neighboring variants as done in Wong et al.

We simply keep these annotations and do a post-filtration later.

isort:skip_file
"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

region_filter = ""
if args["bad_region_expressions"]:
    region_filter = " || " + " || ".join(
        map(
            lambda x: "({})".format(x),
            args["bad_region_expressions"],
        )
    ).replace("$sample_index", args["index_library"])

ShellWrapper(snakemake).run(
    r"""
samples="{args[index_library]}"
samples+=",{args[father]}"
samples+=",{args[mother]}"

# Perform Hard-Filtration --------------------------------------------------------------------------

bcftools view \
    -s "$samples" \
    -i 'INFO/DeNovo[*] == "{args[index_library]}"' \
    -O u \
    {snakemake.input.vcf} \
| bcftools view \
    -e '(ClippedStack == 1) || (FILTER == "Besenbacher"){region_filter}' \
    -O z \
    -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

# Compute Summary ----------------------------------------------------------------------------------

echo -e "chrom\tpos\tid\tref\talt\ttype\tneighbour_samples\tde_novo_origin\tgt_index\tgq_index\tad_index\tgt_father\tgq_father\tad_father\tgt_mother\tgq_mother\tad_mother" \
> {snakemake.output.summary}

bcftools query \
    -f "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%Neighbor\t%DeNovoOrigin[\t%GT\t%GQ\t%AD]\n" \
    {snakemake.output.vcf} \
| {{ grep '^[1-9]' || true; }} \
>> {snakemake.output.summary}
"""
)

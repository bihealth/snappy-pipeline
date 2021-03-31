# -*- coding: utf-8 -*-
"""Snakemake wrapper for hard-filtering the de novo results.

Apply all hard-filters with the following exceptions :

- don't remove variants with dbSNP IDs
- don't remove the neighboring variants as done in Wong et al.

We simply keep these annotations and do a post-filtration later.

isort:skip_file
"""

import os
import sys

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snakemake.shell import shell

region_filter = ""
if snakemake.config["step_config"]["variant_denovo_filtration"]["bad_region_expressions"]:
    region_filter = (
        " || "
        + " || ".join(
            map(
                lambda x: "({})".format(x),
                snakemake.config["step_config"]["variant_denovo_filtration"][
                    "bad_region_expressions"
                ],
            )
        ).replace("$sample_index", snakemake.wildcards.index_library)
    )

shell(
    r"""
set -x
set -euo pipefail

samples="{snakemake.wildcards.index_library}"
samples+=",{snakemake.params.args[father]}"
samples+=",{snakemake.params.args[mother]}"

# Perform Hard-Filtration --------------------------------------------------------------------------

bcftools view \
    -s "$samples" \
    -i 'INFO/DeNovo[*] == "{snakemake.wildcards.index_library}"' \
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

# Create Checksum Files ----------------------------------------------------------------------------

pushd $(dirname {snakemake.output.vcf})

md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi}).md5

md5sum $(basename {snakemake.output.summary}) >$(basename {snakemake.output.summary}).md5
"""
)

# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for frequency filter for variant_filtration.
"""

import os
import sys

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Prelude -----------------------------------------------------------------------------------------

shell.executable("/bin/bash")
shell.prefix("set -eu -o pipefail -x; ")

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

# Short-circuit in case of performing no filtration
if snakemake.wildcards.frequency == "freq_all":
    shell(
        r"""
    # Frequency set to "freq_all", just copy out the data.
    cp -L {snakemake.input.vcf} {snakemake.output.vcf}
    cp -L {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
    """
    )
    sys.exit(0)  # everything went well!


# Actual Filtration -------------------------------------------------------------------------------

# Get shortcut to frequencies set.
frequencies = snakemake.config["step_config"]["variant_filtration"]["frequencies"]

shell(
    r"""
set -x

## do we have to separate de novo and dom freq? I don't think so.
if [[ "{snakemake.wildcards.frequency}" == dominant_freq ]]; then
    # Perform filtration for de novo variants
    include='((GNOMAD_GENOMES_AC_ALL <= {frequencies[ac_dominant]}) | (GNOMAD_GENOMES_AC_ALL = "."))'
    include+=' && ((GNOMAD_GENOMES_AF_POPMAX < {frequencies[af_dominant]}) | (GNOMAD_GENOMES_AF_POPMAX = "."))'
    include+=' && ((GNOMAD_EXOMES_AF_ALL < {frequencies[af_dominant]}) | (GNOMAD_EXOMES_AF_ALL = "."))'
    include+=' && ((GNOMAD_EXOMES_AF_POPMAX < {frequencies[af_dominant]}) | (GNOMAD_EXOMES_AF_POPMAX = "."))'
    bcftools view \
        -i "$include" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
elif [[ "{snakemake.wildcards.frequency}" == recessive_freq ]]; then
    # Perform filtration for homozygous variants, heterozygous in both father and mother
    include='((GNOMAD_GENOMES_AF_ALL < {frequencies[af_recessive]}) | (GNOMAD_GENOMES_AF_ALL = "."))'
    include+=' && ((GNOMAD_GENOMES_AF_POPMAX < {frequencies[af_recessive]}) | (GNOMAD_GENOMES_AF_POPMAX = "."))'
    include+=' && ((GNOMAD_EXOMES_AF_ALL < {frequencies[af_recessive]}) | (GNOMAD_EXOMES_AC_ALL = "."))'
    include+=' && ((GNOMAD_EXOMES_AF_POPMAX < {frequencies[af_recessive]}) | (GNOMAD_EXOMES_AF_POPMAX = "."))'
    include+=' && (DBSNP_COMMON != 1)'
    bcftools view \
        -i "$include" \
        -O z \
        -o {snakemake.output.vcf} \
        {snakemake.input.vcf}
else
    links=1
    ln -sr {snakemake.input.vcf} {snakemake.output.vcf}
    ln -sr {snakemake.input.vcf_md5} {snakemake.output.vcf_md5}
    ln -sr {snakemake.input.vcf_tbi} {snakemake.output.vcf_tbi}
    ln -sr {snakemake.input.vcf_tbi_md5} {snakemake.output.vcf_tbi_md5}
fi

if [[ "${{links-0}}" -ne 1 ]]; then
    tabix -f {snakemake.output.vcf}

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
    md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
fi
"""
)

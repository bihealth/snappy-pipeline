# -*- coding: utf-8 -*-
"""Snakemake wrapper for summarising the de novo counts."""

import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.input), file=tmpf)
    tmpf.flush()
    # Actually run the script
    shell(
        r"""
    set -x

    set -euo pipefail

    get_sample()
    {{
        echo $(basename $1 .summary.txt) \
        | sed -e 's/.*de_novos_hard\.//g' \
        | cut -d - -f 1-4
    }}

    get_caller()
    {{
        echo $(basename $1 .summary.txt) \
        | cut -d . -f 2
    }}

    count_snvs()
    {{
        sample=$(get_sample $1)
        awk -v sample=$sample 'BEGIN {{ count = 0; }} ($3 == "." && $6 == "SNP" && $7 !~ sample)' $1 | wc -l
    }}

    count_snvs_paternal()
    {{
        sample=$(get_sample $1)
        awk -v sample=$sample 'BEGIN {{ count = 0; }} ($3 == "." && $6 == "SNP" && $7 !~ sample && $8 ~ (sample "\\|paternal"))' $1 | wc -l
    }}

    count_snvs_maternal()
    {{
        sample=$(get_sample $1)
        awk -v sample=$sample 'BEGIN {{ count = 0; }} ($3 == "." && $6 == "SNP" && $7 !~ sample && $8 ~ (sample "\\|maternal"))' $1 | wc -l
    }}

    count_indels()
    {{
        sample=$(get_sample $1)
        awk -v sample=$sample 'BEGIN {{ count = 0; }} ($3 == "." && $6 == "INDEL" && $7 !~ sample)' $1 | wc -l
    }}

    (
    echo -e "sample\tcaller\tdenovo_snv\tdenovo_snv_paternal\tdenovo_snv_maternal\tdenovo_indel\tfile"
    for txt in $(sort {tmpf.name}); do
        echo -e "$(get_sample $txt | sed -e 's/-N1-DNA1-WGS1//' -e 's/_/-/g')\t$(get_caller $txt)\t$(count_snvs $txt)\t$(count_snvs_paternal $txt)\t$(count_snvs_maternal $txt)\t$(count_indels $txt)\t$txt"
    done
    ) \
    > {snakemake.output.txt}

    pushd $(dirname {snakemake.output.txt})
    md5sum $(basename {snakemake.output.txt}) >$(basename {snakemake.output.txt}).md5
    """
    )

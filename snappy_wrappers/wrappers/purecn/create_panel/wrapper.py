# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing PureCN panel of normals
"""

import tempfile

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

shell.executable("/bin/bash")

with tempfile.NamedTemporaryFile("wt") as tmpf:
    print("\n".join(snakemake.input.normals), file=tmpf)
    tmpf.flush()

    shell(
        r"""
    set -x

    # Also pipe everything to log file
    if [[ -n "{snakemake.log.log}" ]]; then
        if [[ "$(set +e; tty; set -e)" != "" ]]; then
            rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
            exec &> >(tee -a "{snakemake.log.log}" >&2)
        else
            rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
            echo "No tty, logging disabled" >"{snakemake.log.log}"
        fi
    fi

    # Write out information about conda installation.
    conda list >{snakemake.log.conda_list}
    conda info >{snakemake.log.conda_info}
    md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
    md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

    # Setup auto-cleaned tmpdir
    export tmpdir=$(mktemp -d)
    trap "rm -rf $tmpdir" EXIT

    export R_LIBS_USER=$(dirname {snakemake.input.packages})

    # Find PureCN scripts in extdata
    pkg_folders=$(R --quiet --vanilla -e 'cat(.libPaths(), sep="\n")' | grep -v '^>')
    PURECN=$(for folder in $pkg_folders ; do ls -1 $folder | grep -E '^PureCN$' | sed -e "s#^#$folder/#" ; done)
    PURECN="$PURECN/extdata"

    # Create panel
    Rscript $PURECN/NormalDB.R \
        --out-dir $(dirname {snakemake.output.rds}) \
        --coverage-files {tmpf.name} \
        --genome {config[genome_name]} --assay {config[enrichment_kit_name]}

    # Rename output
    d=$(dirname {snakemake.output.rds})
    mv $d/normalDB_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.rds}
    mv $d/low_coverage_targets_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.bed}
    mv $d/interval_weights_{config[enrichment_kit_name]}_{config[genome_name]}.png {snakemake.output.png}

    # MD5 checksum for main result only
    pushd $(dirname {snakemake.output.rds})
    f=$(basename {snakemake.output.rds})
    md5sum $f > $f.md5
    popd
    """
    )

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

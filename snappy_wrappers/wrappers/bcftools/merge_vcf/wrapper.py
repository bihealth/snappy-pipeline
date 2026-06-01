# -*- coding: utf-8 -*-
"""Wrapper for running bcftools merge - Structural VCF files (CNV, SV)."""

import tempfile
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})
merge_option = args["merge_option"]
gvcf_option = args["gvcf_option"]
sample_names = args["sample_names"]
input_ = args["input"]

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(input_), file=tmpf)
    tmpf.flush()
    ShellWrapper(snakemake).run(
        r"""
    # Method checks if VCF contains sample
    check_vcf() {{
        # Variables
        # vcf=$1
        # sample=$2

        # Check
        if bcftools query --list-samples $1 | grep --quiet --word-regexp $2; then
           return 0
        else
            echo "VCF header doesn't contain sample '$2': $1"
            echo "Samples:" $(bcftools query --list-samples $1)
            exit 1
        fi
    }}

    # Symlink input VCF files
    i=0
    mkdir $TMPDIR/cwd
    for x in $(cat {tmpf.name}); do
        let "i=$i+1"
        ln -s $(readlink -f $x) $TMPDIR/cwd/$i.vcf.gz
        ln -s $(readlink -f $x).tbi $TMPDIR/cwd/$i.vcf.gz.tbi
    done

    # -----------
    # Merge VCFs
    # -----------

    # Define merge option
    merge_option="--merge none"
    if [[ "{merge_option}" != "None" ]]; then
        merge_option="--merge {merge_option}"
    fi

    # Set merge gVCF option
    gvcf_option=""
    if [[ "{gvcf_option}" != "False" ]]; then
        gvcf_option="--gvcf"
    fi

    # If a single sample, there is no need to merge.
    # ``$i`` is reused from previous VCFs to temp dir for-loop.
    if [[ $i -eq 1 ]]; then
        # Validate VCF: contains all expected samples
        while read sample; do
            check_vcf $TMPDIR/cwd/1.vcf.gz $sample
        done < <(echo {sample_names})
        # Copy
        cp $TMPDIR/cwd/1.vcf.gz {snakemake.output.vcf}
        cp $TMPDIR/cwd/1.vcf.gz.tbi {snakemake.output.vcf_tbi}
    else
        out=$(realpath {snakemake.output.vcf})
        pushd $TMPDIR/cwd
        bcftools merge \
            $merge_option $gvcf_option \
            --missing-to-ref \
            --output-type z \
            --output $out \
            *.vcf.gz
        popd
        tabix -f {snakemake.output.vcf}
    fi

    """
    )

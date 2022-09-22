# -*- coding: utf-8 -*-
"""Wrapper for running bcftools merge - VCF files.
"""

import tempfile

from snakemake.shell import shell

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.params.args["input"]), file=tmpf)
    tmpf.flush()
    # Actually run the script.
    shell(
        r"""
    # -----------------------------------------------------------------------------
    # Redirect stderr to log file by default and enable printing executed commands
    exec &> >(tee -a "{snakemake.log.log}")
    set -x
    # -----------------------------------------------------------------------------
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    # Write out information about conda installation
    conda list > {snakemake.log.conda_list}
    conda info > {snakemake.log.conda_info}

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
    if [[ "{snakemake.params.args[merge_option]}" != "None" ]]; then
        merge_option="--merge {snakemake.params.args[merge_option]}"
    fi

    # Set merge gVCF option
    gvcf_option=""
    if [[ "{snakemake.params.args[gvcf_option]}" != "False" ]]; then
        gvcf_option="--gvcf"
    fi

    # If a single sample, there is no need to merge.
    # ``$i`` is reused from previous BCFs for-loop.
    if [[ $i -eq 1 ]]; then
        # Validate VCF: contains all expected samples
        while read sample; do
            check_vcf $TMPDIR/cwd/1.vcf.gz $sample
        done < <(echo {snakemake.params.args[sample_names]})
        # Copy
        cp $TMPDIR/cwd/1.vcf.gz {snakemake.output.vcf}
        cp $TMPDIR/cwd/1.vcf.gz.tbi {snakemake.output.tbi}
    else
        out=$(realpath {snakemake.output.vcf})
        pushd $TMPDIR/cwd
        bcftools merge \
            $merge_option $gvcf_option \
            --output-type z \
            --output $out \
            *.vcf.gz
        popd
        tabix -f {snakemake.output.vcf}
    fi

    pushd $(dirname {snakemake.output.vcf})
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf_md5})
    md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi_md5})
    """
    )

# Compute MD5 sums of logs
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
"""
)

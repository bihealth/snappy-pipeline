# -*- coding: utf-8 -*-
"""Wrapper for running Platypus in joint mode."""

import os

from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

args = getattr(snakemake.params, "args", {})
args_ignore_chroms = ""
if ignore_chroms := args.get("ignore_chroms"):
    args_ignore_chroms = " ".join(["--ignore-chroms"] + ignore_chroms)

reference_path = snakemake.input.reference
num_threads = args["num_threads"]
split_complex_mnvs = str(args["split_complex_mnvs"])


split_script = os.path.join(os.path.dirname(__file__), "splitMNPsAndComplex.py")

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

# Create temporary directory and organize cleanup
export TMPDIR=$(mktemp -d)
trap "rm -rf \"$TMPDIR\"" EXIT

# Platypus Variant Calling ----------------------------------------------------

REF={reference_path}
out_final={snakemake.output.vcf}
#out_tmp=$TMPDIR/out_tmp.vcf.gz
out_tmp=${{out_final%.vcf.gz}}.tmp.vcf.gz

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

platypus callVariants \
    --logFileName=$(dirname {snakemake.log})/platypus.log \
    --bamFiles=$(echo "{snakemake.input.bam}" | tr ' ' ',') \
    --nCPU={num_threads} \
    --refFile=$REF \
    --output=${{out_tmp%.gz}} \
    --regions=$(snappy-genome_windows \
                    --fai-file $REF.fai \
                    --window-size 1000000000 \
                    {args_ignore_chroms} \
                | sed -e 's/,//g' \
                | tr '\n' ',' \
                | sed -e 's/,$//g')

if [[ "{split_complex_mnvs}" == "True" ]]; then
    cat ${{out_tmp%.gz}} \
    | python2 {split_script} \
    | snappy-vcf_sort $REF.fai \
    > ${{out_final%.gz}}
else
    mv ${{out_tmp%.gz}} ${{out_final%.gz}}
fi

bgzip ${{out_final%.gz}}
tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) > $(basename {snakemake.output.vcf_tbi}).md5
"""
)

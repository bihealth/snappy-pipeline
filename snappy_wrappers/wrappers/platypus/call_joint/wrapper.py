# -*- coding: utf-8 -*-
"""Wrapper for running Platypus in joint mode.
"""

import os

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

args_ignore_chroms = ""
if snakemake.params.args["ignore_chroms"]:
    args_ignore_chroms = " ".join(["--ignore-chroms"] + snakemake.params.args["ignore_chroms"])

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

REF={snakemake.config[static_data_config][reference][path]}
out_final={snakemake.output.vcf}
#out_tmp=$TMPDIR/out_tmp.vcf.gz
out_tmp=${{out_final%.vcf.gz}}.tmp.vcf.gz

# Hack: get back bin directory of base/root environment.
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

platypus callVariants \
    --logFileName=$(dirname {snakemake.log})/platypus.log \
    --bamFiles=$(echo "{snakemake.input.bam}" | tr ' ' ',') \
    --nCPU={snakemake.config[step_config][somatic_variant_calling][platypus_joint][num_threads]} \
    --refFile=$REF \
    --output=${{out_tmp%.gz}} \
    --regions=$(snappy-genome_windows \
                    --fai-file $REF.fai \
                    --window-size 1000000000 \
                    {args_ignore_chroms} \
                | sed -e 's/,//g' \
                | tr '\n' ',' \
                | sed -e 's/,$//g')

if [[ "{snakemake.config[step_config][somatic_variant_calling][platypus_joint][split_complex_mnvs]}" == "True" ]]; then
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
md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi}).md5
"""
)

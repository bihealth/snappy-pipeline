# -*- coding: utf-8 -*-
"""Wrapper for running SV2 re-genotyping step on merged ERDS calls.
"""

import os.path

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bih-charite.de"

# Get path to this file's (wrapper.py) directory.
base_dir = os.path.dirname(os.path.realpath(__file__))

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
# exec &> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

# We need to use the temporary directory on the local disk here as a bazillion
# files are created by SV2.
export TMPDIR=/tmp

export LC_ALL=C
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

mkdir -p $TMPDIR/{{tmp,out}}

vcf=$(basename {snakemake.output.vcf} .gz)

# TODO: put correct sex here
echo -e "FAM\t{snakemake.wildcards.library_name}\t0\t0\t2\t2" \
>$TMPDIR/tmp/tmp.ped

config_ref={snakemake.config[static_data_config][reference][path]}
if [[ "$config_ref" == *h37* ]] || [[ "$config_ref" == *hs37* ]] || [[ "$config_ref" == *hg19* ]]; then
    arg_ref=hg19
elif [[ "$config_ref" == *h38* ]] || [[ "$config_ref" == *h?38* ]]; then
    arg_ref=hg38
elif [[ "$config_ref" == *m38* ]] || [[ "$config_ref" == *mm10* ]]; then
    arg_ref=mm10
else
    >&2 echo 'Cannot guess reference name. Giving up!'
    exit 1
fi

cat >$TMPDIR/sv2.ini <<"EOF"
[FASTA_PATHS]
hg19 = {snakemake.config[step_config][wgs_cnv_calling][sv2][path_hg19]}
hg38 = {snakemake.config[step_config][wgs_cnv_calling][sv2][path_hg38]}
mm10 = {snakemake.config[step_config][wgs_cnv_calling][sv2][path_mm10]}

[INSTALL_DIR]
sv2_home = None

[RESOURCE_DIR]
sv2_resource = {snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_resource]}
EOF

cat >$TMPDIR/sv2_clf.json <<"EOF"
{{
    "default": {{
        "delmsc": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_highcov_del_malesexchrom_svm_clf.pkl",
        "duphar": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_highcov_dup_har_svm_clf.pkl",
        "dellt": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_highcov_del_lt1kb_svm_clf.pkl",
        "dupbrk": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_lowcov_dup_breakpoint_svm_clf.pkl",
        "dupmsc": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_lowcov_dup_malesexchrom_svm_clf.pkl",
        "delgt": "{snakemake.config[step_config][wgs_cnv_calling][sv2][path_sv2_models]}/1kgp_highcov_del_gt1kb_svm_clf.pkl"
    }}
}}
EOF

head -n 1000 $TMPDIR/sv2.ini $TMPDIR/sv2_clf.json

sv2 \
    -ini $TMPDIR/sv2.ini \
    -bam {snakemake.input.bam} \
    -vcf {snakemake.input.vcf_cnv} \
    -snv {snakemake.input.vcf_small} \
    -ped $TMPDIR/tmp/tmp.ped \
    -g $arg_ref \
    -o $vcf \
    -O $TMPDIR/out \
    -tmp-dir $TMPDIR/tmp

sed \
    -e 's/ID=GL,Number=2|3/ID=GL,Number=G/g' \
    -e 's/GAP>,/GAP,/g' \
    $TMPDIR/out/sv2_genotypes/$vcf \
| awk \
    -F $'\t' \
    'BEGIN {{ OFS=FS }}
     /^#/ {{ print }}
     (!/^#/) {{ $4 = 'N'; gsub(",", ";", $7); print }}' \
| bcftools annotate \
    --set-id +'%CHROM\_%POS\_%INFO/END\_%SVTYPE' \
    --remove "FILTER/GENOTYPEFAIL,FILTER/NOALT" \
    -O z \
    -o {snakemake.output.vcf}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf}).tbi >$(basename {snakemake.output.vcf}).tbi.md5
"""
)

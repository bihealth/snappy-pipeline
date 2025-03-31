# -*- coding: utf-8 -*-
"""Wrapper code for GATK PhaseByTransmission"""

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.prefix("set -euo pipefail; ")

args = getattr(snakemake.params, "args", {})

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

which tabix
which bcftools
which gatk_nonfree

# Extract pedigree only
members=$(cut -f 2 {snakemake.input.ped} | tr '\n' ',' | sed -e 's/,$//'g)
bcftools view \
    --threads 4 \
    -s "$members" \
    {snakemake.input.vcf} \
    -O z \
    -o $TMPDIR/trio_only.vcf.gz
tabix -f $TMPDIR/trio_only.vcf.gz

# Call wrapper
MALLOC_ARENA_MAX=4 \
gatk_nonfree \
    -Xmx10g \
    -Djava.io.tmpdir=$TMPDIR \
    --analysis_type PhaseByTransmission \
    -nct 1 \
    --pedigreeValidationType SILENT \
    --FatherAlleleFirst \
    --DeNovoPrior {args[de_novo_prior]} \
    --pedigree {snakemake.input.ped} \
    --variant $TMPDIR/trio_only.vcf.gz \
    --out {snakemake.output.vcf} \
    --reference_sequence {args[reference]}

tabix -f {snakemake.output.vcf}

pushd $(dirname {snakemake.output.vcf})
md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf}).md5
md5sum $(basename {snakemake.output.vcf_tbi}) >$(basename {snakemake.output.vcf_tbi}).md5
popd
"""
)

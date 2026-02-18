# -*- coding: utf-8 -*-
"""Wrapper for creating personalized proteome from a set of (germline) variants"""

import os

from snakemake.shell import shell

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

script = os.path.join(os.path.dirname(__file__), "create_proteome.R")

add_reference = "reference" if snakemake.input.get("add_reference", False) else "NULL"
path_proteome = snakemake.input.get("path_proteome", "")


shell(
    r"""
set -x

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

tmp=$(mktemp -d)
vcf=$tmp/normalized.vcf.gz

bcftools norm \
    --fasta-file {snakemake.input.reference} \
    --remove-duplicates --multiallelics -both \
    --output-type z --output $vcf --write-index=tbi \
    {snakemake.input.vcf}

cat << __EOF > $tmp/script.R
source("{script}")

tx <- GenomicFeatures::makeTxDbFromGFF("{snakemake.input.features}")
reference <- Biostrings::readDNAStringSet("{snakemake.input.reference}")
variants <- "$vcf"
stopifnot(all(lengths(VariantAnnotation::alt(variants)) == 1))

reference <- clean_genome(reference, GenomeInfoDb::seqinfo(tx))

genomes <- insert_variants(reference, variants)

proteins <- donor_proteome(genomes, tx, reference={add_reference})

if ("{path_proteome}" != "") proteins <- unique(c(proteins, as.character(Biostrings::readDNAStringSet("{path_proteome}"))))

Biostrings::writeXStringSet(proteins, "{snakemake.output.proteins}")
__EOF

R --vanilla -e $tmp/script.R

pushd $(dirname {snakemake.output.proteins})
f=$(basename {snakemake.output.proteins})
md5sum $f > $f.md5
popd
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

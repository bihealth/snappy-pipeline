import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Optionally get path to coverage VCF file.
coverage_vcf = " ".join(getattr(snakemake.input, "vcf_cov", []))

args = getattr(snakemake.params, "args", {})

# Get shortcut to "fix_manta_invs.py" postprocessing script
fix_manta_invs = os.path.join(
    os.path.dirname(__file__),
    "fix_manta_invs.py",
)

ShellWrapper(snakemake).run(
    r"""
set -x

samples=$(cut -f 2 {snakemake.input.ped} | tr '\n' ',' | sed -e 's/,$//g')

# Fix the Manta inversions
i=0
for vcf in {snakemake.input.vcf}; do
    let "i=$i+1"
    num=$(printf %03d $i)

    python3 {fix_manta_invs} \
        --reference-fasta {args[reference]} \
        --input-vcf $vcf \
        --output-vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf
    bcftools sort -o $TMPDIR/fixed_bnd_to_inv.$num.vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf

    # Fixup SVLEN=1 to SVLEN=.
    sed -i -e 's/ID=SVLEN,Number=1/ID=SVLEN,Number=./g' $TMPDIR/fixed_bnd_to_inv.$num.vcf
    # Fixup MELT header
    sed -i -e "s/seperated by '..'/separated by '\\\\\\\\|'/" $TMPDIR/fixed_bnd_to_inv.$num.vcf

    # Add the missing "GT" tag
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    > $TMPDIR/header.gt.txt

    # Annotate and normalise the VCF files.
    bcftools annotate \
        -h $TMPDIR/header.gt.txt \
        $TMPDIR/fixed_bnd_to_inv.$num.vcf |
        bcftools norm \
            -m -any \
            -c w \
            --fasta-ref {args[reference]} \
            -O z \
            -o $TMPDIR/final_for_import.$num.vcf.gz
    tabix -s1 -b2 -e2 -f $TMPDIR/final_for_import.$num.vcf.gz
done

cat <<"EOF" > $TMPDIR/feature-effects.tsv
case_id
set_id
sv_uuid
refseq_gene_id
refseq_transcript_id
refseq_transcript_coding
refseq_effect
ensembl_gene_id
ensembl_transcript_id
ensembl_transcript_coding
ensembl_effect
EOF

cat $TMPDIR/feature-effects.tsv \
| tr '\n' '\t' \
| sed -e 's/\t$/\n/g' \
| gzip \
>{snakemake.output.feature_effects}

# Perform Mehari structural variant annotation.
mehari \
    annotate \
    strucvars \
    --path-input-ped {snakemake.input.ped} \
    $(for p in $TMPDIR/final_for_import.*.vcf.gz; do \
        echo --path-input-vcf $p; \
    done) \
    --path-output-tsv >(gzip -c > {snakemake.output.gts})

cat <<EOF | gzip -c > {snakemake.output.db_infos}
genomebuild	db_name	release
GRCh37	clinvar	20210728
GRCh37	gnomad_exomes	r2.1.1
GRCh37	gnomad_genomes	r2.1.1
EOF
"""
)

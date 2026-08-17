import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

# Optionally get path to coverage VCF file.
coverage_vcf = " ".join(getattr(snakemake.input, "vcf_cov", []))

args = getattr(snakemake.params, "args", {})
export_config = args["config"]

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
        --reference-fasta {snakemake.input.reference} \
        --input-vcf $vcf \
        --output-vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf
    bcftools sort -o $TMPDIR/fixed_bnd_to_inv.$num.vcf $TMPDIR/fixed_bnd_to_inv_unsorted.$num.vcf

    # Add the missing "GT" tag
    echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' \
    > $TMPDIR/header.gt.txt

    bcftools annotate \
        -h $TMPDIR/header.gt.txt \
        $TMPDIR/fixed_bnd_to_inv.$num.vcf \
        -O z \
        -o $TMPDIR/final_for_import.$num.vcf.gz
    tabix -s1 -b2 -e2 -f $TMPDIR/final_for_import.$num.vcf.gz
done

# Compatibility mode with currently deployed VarFish Server
compatibility_option="--opt-out callers-array"

# Execute VarFish Annotator
varfish-annotator \
    annotate-svs \
    -XX:MaxHeapSize=10g \
    -XX:+UseG1GC \
    \
    --release {export_config[release]} \
    \
    $(if [[ "{coverage_vcf}" != "" ]]; then \
        for path in {coverage_vcf}; do \
            echo --coverage-vcf $path; \
        done; \
    fi) \
    \
    --db-path {export_config[path_db]} \
    --refseq-ser-path {export_config[path_refseq_ser]} \
    --ensembl-ser-path {export_config[path_ensembl_ser]} \
    --input-ped {snakemake.input.ped} \
    \
    $(for vcf in $TMPDIR/final_for_import.*.vcf.gz; do \
        echo --input-vcf $vcf; \
    done) \
    --output-db-info {snakemake.output.db_infos} \
    --output-gts {snakemake.output.gts} \
    --output-feature-effects {snakemake.output.feature_effects} \
    $compatibility_option
"""
)

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
export_config = args["config"]

ShellWrapper(snakemake).run(
    r"""
set -x

# Extract around BED file, if given.
if [[ -n "{export_config[path_exon_bed]}" ]] && [[ "{export_config[path_exon_bed]}" != "None" ]]; then
    set -e
    bcftools view \
        -R {export_config[path_exon_bed]} \
        {snakemake.input.vcf} \
    | bcftools sort -T $TMPDIR \
    | bcftools norm -d all \
    | bgzip -c \
    > $TMPDIR/tmp.vcf.gz
    tabix -f $TMPDIR/tmp.vcf.gz
else
    set -e
    ln -sr {snakemake.input.vcf} $TMPDIR/tmp.vcf.gz
    ln -sr {snakemake.input.vcf}.tbi $TMPDIR/tmp.vcf.gz.tbi
fi

# Execute VarFish Annotator
varfish-annotator \
    annotate \
    -XX:MaxHeapSize=10g \
    -XX:+UseG1GC \
    \
    --release {export_config[release]} \
    \
    --self-test-chr1-only \
    --ref-path {snakemake.input.reference} \
    --db-path {export_config[path_db]} \
    --refseq-ser-path {export_config[path_refseq_ser]} \
    --ensembl-ser-path {export_config[path_ensembl_ser]} \
    --input-ped {snakemake.input.ped} \
    \
    --input-vcf $TMPDIR/tmp.vcf.gz \
    --output-db-info {snakemake.output.db_infos} \
    --output-gts {snakemake.output.gts}

# Copy out PED file to output
cp -H {snakemake.input.ped} {snakemake.output.ped}
"""
)

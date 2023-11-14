from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

shell.executable("/bin/bash")

# Build standard filters argument
arg_std_filters = ""
if snakemake.config["step_config"]["variant_calling"]["freebayes"]["use_standard_filters"]:
    arg_std_filters = "--standard-filters"

# Build ignore regions argument
ignore_chroms = snakemake.config["step_config"]["variant_calling"]["freebayes"].get(
    "ignore_chroms", ""
)
if ignore_chroms:
    arg_ignore_chroms = "--ignore-chroms " + " ".join(map(repr, ignore_chroms))
else:
    arg_ignore_chroms = ""

# Get number of threads parameter
num_threads = snakemake.config["step_config"]["variant_calling"]["freebayes"]["num_threads"]

# Define which freebayes mode to use
freebayes = "freebayes"
if num_threads > 1:
    freebayes_parallel_tpl = (
        "freebayes-parallel <(fasta_generate_regions.py {ref} {window_length}) {threads} "
    )
    freebayes = freebayes_parallel_tpl.format(
        ref=snakemake.config["static_data_config"]["reference"]["path"],
        window_length=snakemake.params.window_length,
        threads=num_threads,
    )

shell(
    r"""
set -x
# Write out information about conda installation.
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty,logging disabled" >"{snakemake.log.log}"
    fi
fi
# Export library
export LD_LIBRARY_PATH=$(dirname $(which bgzip))/../lib
# Hack: get back bin directory of base/root environment
export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
# Create shortcut to reference
export REF={snakemake.config[static_data_config][reference][path]}
# Call Freebayes
{freebayes} \
    --fasta-reference $REF \
    --genotype-qualities \
    --min-repeat-entropy {snakemake.params.min_repeat_entropy} \
    --haplotype-length {snakemake.params.haplotype_length} \
    --min-alternate-fraction {snakemake.params.min_alternate_fraction} \
    --min-mapping-quality {snakemake.params.min_mapping_quality} \
    $(echo {snakemake.input} | tr ' ' '\n' | grep '\.bam$') \
| bcftools norm --fasta-ref $REF \
| snappy-vcf_sort $REF.fai \
> $TMPDIR/tmp.vcf
# Filter VCF if requested
if [[ -n "{arg_ignore_chroms}" ]]; then
    # Create include regions bed file
    snappy-genome_windows --fai-file $REF.fai {arg_ignore_chroms} \
        | sed "s/,//g" | sed "s/[:-]/\t/g" \
        > $TMPDIR/tmp.bed
    # Intersect VCF file
    bedtools intersect -a  $TMPDIR/tmp.vcf -b $TMPDIR/tmp.bed -wa -header \
        | bgzip -c \
        > {snakemake.output.vcf}
else
    bgzip -c $TMPDIR/tmp.vcf > {snakemake.output.vcf}
fi
# compute tabix index of the resulting VCF file
tabix -f {snakemake.output.vcf}
pushd $(dirname {snakemake.output.vcf}) && \
    md5sum $(basename {snakemake.output.vcf}) > $(basename {snakemake.output.vcf}).md5 && \
    md5sum $(basename {snakemake.output.tbi}) > $(basename {snakemake.output.tbi}).md5
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} > {snakemake.log.log_md5}
"""
)

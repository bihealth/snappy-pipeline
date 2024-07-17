import os

from snakemake.shell import shell

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]

algorithms = snakemake.params.args["algorithms"]
epitope_length = snakemake.params.args["epitope_length"]
normal_smaple = snakemake.params.args["normal_sample"]
BINDING_THRESHOLD = f"-b {snakemake.params.args["BINDING_THRESHOLD"]} "
percentile_threshold = (
    f"--percentile-threshold {snakemake.params.args["percentile_threshold"]} "
    if snakemake.params.args["percentile_threshold"] != None
    else ""
)


allele_specific_binding_thresholds = (
    f"--allele-specific-binding-thresholds "
    if snakemake.params.args["allele_specific_binding_thresholds"]
    else ""
)

aggregate_inclusion_binding_threshold = (
    f"--aggregate-inclusion-binding-threshold {snakemake.params.args["aggregate_inclusion_binding_threshold"]} ",
)

netmhc_stab = f"--netmhc-stab " if snakemake.params.args["netmhc_stab"] else ""

NET_CHOP_THRESHOLD = (f"--net-chop-threshold {snakemake.params.args["NET_CHOP_THRESHOLD"]} ",)

PROBLEMATIC_AMINO_ACIDS = (
    f"--problematic-amino-acids {snakemake.params.args["PROBLEMATIC_AMINO_ACIDS"]} ",
)

# if {snakemake.params.args["run_reference_proteome_similarity"]}:
#     run_reference_proteome_similarity= f"--run-reference-proteome-similarity",

FAST_SIZE = (f"-s {snakemake.params.args["FAST_SIZE"]} ",)

exclude_NAs = "--exclude-NAs " if snakemake.params.args["exclude_NAs"] else ""

NORMAL_COV = (f"--normal-cov {snakemake.params.args["NORMAL_COV"]} ",)
TDNA_COV = (f"--tdna-cov {snakemake.params.args["TDNA_COV"] }",)
TRNA_COV = (f"--trna-cov {snakemake.params.args["TRNA_COV"]} ",)
NORMAL_VAF = (f"--normal-vaf {snakemake.params.args["NORMAL_VAF"]} ",)

maximum_transcript_support_level = (
    f"--maximum-transcript-support-level {snakemake.params.args["maximum_transcript_support_level"]} ",
)

pass_only = "--pass-only " if snakemake.params.args["pass_only"] else ""

tumor_purity = (
    f"--tumor-purity {snakemake.params.args["tumor_purity"]} "
    if snakemake.params.args["tumor_purity"] != None
    else ""
)

op_dir = "/".join(snakemake.output.all_epitopes.split("/")[:-2])

files_to_bind = {
    "combine_vcf": snakemake.input.combine_vcf,
    "op_dir": op_dir,
}

# TODO: Put the following in a function (decide where...)
# Replace with full absolute paths
files_to_bind = {k: os.path.realpath(v) for k, v in files_to_bind.items()}
# Directories that mut be bound
dirs_to_bind = {k: os.path.dirname(v) for k, v in files_to_bind.items()}
# List of unique directories to bind: on cluster: <directory> -> from container: /bindings/d<i>)
bound_dirs = {e[1]: e[0] for e in enumerate(list(set(dirs_to_bind.values())))}
# Binding command
bindings = " ".join(["-B {}:/bindings/d{}".format(k, v) for k, v in bound_dirs.items()])
bound_files = {
    k: "/bindings/d{}/{}".format(bound_dirs[dirs_to_bind[k]], os.path.basename(v))
    for k, v in files_to_bind.items()
}

shell.executable("/bin/bash")

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
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
# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Compute md5 checksum
md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo "$checksum"
}}

set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation
conda list > {snakemake.log.conda_list}
conda info > {snakemake.log.conda_info}

hla_types=$(cat {snakemake.input.hla_normal_dna} {snakemake.input.hla_tumor_dna} {snakemake.input.hla_tumor_rna} | sed 's/^/HLA-/' | sort | uniq | tr '\n' ',' | sed 's/,$//')
cmd="pvacseq run --normal-sample-name {normal_smaple} \
    -e1 {epitope_length} \
    --iedb-install-directory /opt/iedb \
    -t {snakemake.threads} \
    {bound_files[combine_vcf]} \
    {snakemake.wildcards.tumor_library} \
    $hla_types \
    {algorithms} \
    {bound_files[op_dir]} \
    {BINDING_THRESHOLD}\
    {percentile_threshold}\
    {allele_specific_binding_thresholds}\
    {aggregate_inclusion_binding_threshold}\
    {NET_CHOP_THRESHOLD}\
    {PROBLEMATIC_AMINO_ACIDS}\
    {FAST_SIZE}\
    {NORMAL_COV}\
    {TDNA_COV}\
    {TRNA_COV}\
    {NORMAL_VAF}\
    {exclude_NAs}\
    {maximum_transcript_support_level}"
echo 'TMPDIR=/bindings/d2' > $TMPDIR/{snakemake.wildcards.tumor_library}.sh
echo $cmd >> $TMPDIR/{snakemake.wildcards.tumor_library}.sh
apptainer exec --home $PWD -B $TMPDIR:/bindings/d2 {bindings} {config[path_container]} bash /bindings/d2/{snakemake.wildcards.tumor_library}.sh

md5 {snakemake.output.all_epitopes} > {snakemake.output.all_epitopes_md5}
md5 {snakemake.output.filtered_epitopes} > {snakemake.output.filtered_epitopes_md5}
"""
)

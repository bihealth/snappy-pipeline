import os

from snakemake import shell
from snakemake.io import Namedlist

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

# TODO: Put the following in a function (decide where...)
def bindings_for_container(files_to_bind: Namedlist, base_container_path: str = "inputs", access: str = "ro") -> tuple[str, dict[str, str]]:
    fns = {}
    for key, path in files_to_bind:
        assert not isinstance(path, list), f"Input files cannot be lists (offender is '{key}')"
        assert key not in fns, f"Duplicate input '{key}'"
        fns[key] = os.path.realpath(path)

    ds = {}
    for key, path in fns.items():
        ds[key] = os.path.dirname(path)

    bindings = {e[1]: e[0] for e in enumerate(list(set(ds.values())))}

    container_fns = {}
    for key, path in fns.items():
        container_fns[key] = f"/{base_container_path}/d{bindings[ds[key]]}/{os.path.basename(path)}"

    binding_cmd = " ".join([f"--bind {d}:/{base_container_path}/d{i}:{access}" for d, i in bindings.items()])

    return (binding_cmd, container_fns)

(input_bindings, input_fns) = bindings_for_container(snakemake.input)
(output_bindings, output_fns) = bindings_for_container(snakemake.output, base_container_path="outputs", access="rw")

if isinstance(args["algorithms"], list):
    algorithms = ",".join(args["algorithms"])
else:
    algorithms = args["algorithms"]
alleles = ",".join(args["alleles"])

pVACseq_args = [f"--netmhciipan-version {args['netmhciipan_version']}"]

if args["class_i_epitope_length"]:
    pVACseq_args.append("--class-i-epitope-length" + ",".join(map(str, args["class_i_epitope_length"])))
if args["class_ii_epitope_length"]:
    pVACseq_args.append("--class-ii-epitope-length" + ",".join(map(str, args["class_ii_epitope_length"])))

if args["percentile_threshold"] is not None:
    pVACseq_args.append(f"--percentile-threshold {args['percentage_threshold']} --percentile-threshold-strategy {args['percentile_threshold_strategy']}")
else:
    pVACseq_args.append(f"--binding-threshold {args['binding_threshold']}")
if args["allele_specific_binding_thresholds"]:
    pVACseq_args.append("--allele-specific-binding-thresholds")

pVACseq_args.append(f"--top-score-metric {args['top_score_metric']}")
pVACseq_args.append(f"--top-score-metric2 {args['top_score_metric2']}")

for key in ("normal_cov", "tdna_cov", "trna_cov", "normal_vaf", "tdna_vaf", "trna_vaf", "minimum_fold_change", "expn_val"):
    pVACseq_key = key.replace("_", "-")
    pVACseq_args.append(f"--{pVACseq_key} {args[key]}")

pVACseq_args.append(f"--transcript-prioritization-strategy {args['transcript_prioritization_strategy']}")
if args["maximum_transcript_support_level"] is not None:
    pVACseq_args.append(f"--maximum-transcript-support-level {args['maximum_transcript_support_level']}")
if args["biotypes"]:
    pVACseq_args.append("--biotypes " + ",".join(args["biotypes"]))
if args["allow_incomplete_transcript"]:
    pVACseq_args.append("--allow-incomplete-transcript")

if args["net_chop_method"] is not None:
    pVACseq_args.append(f"--net-chop-method {args['net_chop_method']}")
if args["netmhc_stab"]:
    pVACseq_args.append("--netmhc-stab")
if args["allele_specific_anchors"]:
    pVACseq_args.append(f"--allele_specific_anchors --anchor-contribution-threshold {args['anchor_contribution_threshold']}")

if args["problematic_amino_acids"]:
    pVACseq_args.append("--problematic-amino-acids" + ",".join(args["problematic_amino_acids"]))

if args["exclude_NAs"]:
    pVACseq_args.append("--exclude-NAs")

for key in ("downstream_sequence_length", "aggregate_inclusion_binding_threshold", "aggregate_inclusion_count_limit"):
    pVACseq_key = key.replace("_", "-")
    pVACseq_args.append(f"--{pVACseq_key} {args[key]}")

if args["run_reference_proteome_similarity"]:
    pVACseq_args.append(f"--run-reference-proteome-similarity --peptide-fasta {input_fns['peptides']}")

if args["genes_of_interest_file"]:
    pVACseq_args.append(f"--genes-of-interest-file {input_fns['genes']}")

pVACseq_args = " \\\n    ".join(pVACseq_args)

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

# Re-home pVACtools to avoid conflicts with user's setting (.bashrc, ...)
# Create home on TMPDIR? But then the calling script is gone...
home=$(dirname {snakemake.output})/../home
home=$(realpath $home)
rm -rf $home
mkdir -p $home

cat << __EOF > $home/run_pVACseq.sh
pvacseq run --n-threads {snakemake.threads} \\
    --normal-sample-name {args[normal_sample]} \\
    --iedb-install-directory /opt/iedb \\
    {pVACseq_args} \\
    {input_fns[vcf]} {args[tumor_sample]} {alleles} \\
    $(dirname {output_fns[done]})
__EOF
chmod +x $home/run_pVACseq.sh

apptainer exec \
    --home $home --bind $TMPDIR:$TMPDIR:rw \
    {input_bindings} {output_bindings} {snakemake.input[container]} bash run_pVACseq.sh
"""
)

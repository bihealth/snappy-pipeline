import csv
import os

from snakemake import shell
from snakemake.io import Namedlist

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})
extra_args = " ".join(args.get("extra_args", []))

# TODO: Extend to lists, check for directories in directories [outputs/rw precedence on inputs/ro] & put in a function (decide where...)
def bindings_for_container(
    files_to_bind: Namedlist,
    base_container_path: str,
    access: str = "ro",
    exclude_bind: list[str] = []
) -> tuple[str, dict[str, str]]:
    fns = {}
    for key, path in files_to_bind.items():
        if key in exclude_bind:
            continue
        assert not isinstance(path, list), f"Input files cannot be lists (offender is '{key}')"
        assert key not in fns, f"Duplicate input '{key}'"
        fns[key] = os.path.realpath(path)

    ds = {}
    for key, path in fns.items():
        ds[key] = os.path.dirname(path)

    bindings = {e[1]: e[0] for e in enumerate(list(set(ds.values())))}

    container_fns = {}
    for key, path in fns.items():
        container_fns[key] = f"{base_container_path}/d{bindings[ds[key]]}/{os.path.basename(path)}"

    binding_cmd = " ".join([f"--bind {d}:{base_container_path}/d{i}:{access}" for d, i in bindings.items()])

    return (binding_cmd, container_fns)

(input_bindings, input_fns) = bindings_for_container(
    snakemake.input, base_container_path="/inputs", exclude_bind=args["exclude_bind"]
)
(output_bindings, output_fns) = bindings_for_container(
    snakemake.output, base_container_path="/outputs", access="rw", exclude_bind=args["exclude_bind"]
)

if peptides := getattr(snakemake.input, "peptides", ""):
    peptides = f"--run-reference-proteome-similarity --peptide-fasta {input_fns['peptides']}"
if genes := getattr(snakemake.input, "genes", ""):
    genes = f"--genes-of-interest-file {input_fns['genes']}"

alleles = []
with open(snakemake.input.alleles, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        alleles.append(row["HLA Allele"])
alleles = ",".join(sorted(list(set(alleles))))

algorithms = " ".join(args["algorithms"])
class_i_epitope_length = ",".join(map(str, args["lengths"]["class_i"]))
class_ii_epitope_length = ",".join(map(str, args["lengths"]["class_ii"]))

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

export LC_ALL=C.UTF-8

# snakemake.output.done is needed because pVACtools modules create MHC_class_iI? folders
# so the path to bind is work/mapper.caller.annotator(.filtered)?.tumor_dna/out
# rather than the dirname of snakemake.output.all (which is the MHC_class_II? subdir of the
# bound directory)
out=$(dirname {snakemake.output.done})
rm -rf $out
mkdir -p $out
out=$(realpath $out)

tmp=$(dirname $out)/tmp/pvacsplice
rm -rf $tmp
mkdir -p $tmp

scripts=$(dirname $out)/scripts
mkdir -p $scripts

cat << __EOF > $scripts/run_pVACsplice.sh
export TMPDIR=/short_tmp
pvacsplice run --n-threads {snakemake.threads} \\
    --normal-sample-name {args[normal_sample]} \\
    --iedb-install-directory /opt/iedb \\
    --class-i-epitope-length {class_i_epitope_length} --class-ii-epitope-length {class_ii_epitope_length} \\
    {extra_args} \\
    {peptides} {genes} \\
    {input_fns[junctions]} \\
    {args[tumor_sample]} {alleles} {algorithms} \\
    $(dirname {output_fns[done]}) \\
    {input_fns[annotated]} {input_fns[reference]} {input_fns[features]}
__EOF
chmod +x $scripts/run_pVACsplice.sh

# pVACsplice (version 7.1.1) fails with an error when there are no predicted neo-epitopes from splice variants
# The current workaround is to prevent snakemake to fail, and to manually create the required directories
# before creating an empty target file.
# Short path tmp to avoid AF_UNIX path too long error
apptainer exec \
    --no-home --bind $tmp:/short_tmp:rw --bind $scripts:/scripts:ro \
    {input_bindings} {output_bindings} {snakemake.input[container]} bash /scripts/run_pVACsplice.sh  \; || true

link_missing() {{
    combined=$1
    cls_i=$2
    cls_ii=$3

    if [[ ! -e $combined ]]
    then
        mkdir -p $(dirname $combined)
        if [[ ! -e $cls_i ]]
        then
            mkdir -p $(dirname $cls_i)
            if [[ ! -e $cls_ii ]]
            then
                mkdir -p $(dirname $cls_ii)
                echo "WARNING: Missing $cls_i and $cls_ii, pVACsplice didn't finish gracefully."
                touch $cls_i
                touch $cls_ii
                touch $combined
                return 0
            fi
            touch $cls_i
            ln -sr $cls_ii $combined
        else
            mkdir -p $(dirname $cls_ii)
            touch $cls_ii
            ln -sr $cls_i $combined
        fi
    fi
}}


link_missing {snakemake.output[all.Combined]} {snakemake.output[all.MHC_I]} {snakemake.output[all.MHC_II]}
link_missing {snakemake.output[aggregated.Combined]} {snakemake.output[aggregated.MHC_I]} {snakemake.output[aggregated.MHC_II]}
link_missing {snakemake.output[filtered.Combined]} {snakemake.output[filtered.MHC_I]} {snakemake.output[filtered.MHC_II]}

ls -altrR $(dirname {snakemake.output.done})
touch {snakemake.output.done}
"""
)

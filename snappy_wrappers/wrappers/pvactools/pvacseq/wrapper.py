import os

from snakemake import shell
from snakemake.io import Namedlist

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

# TODO: Put the following in a function (decide where...)
def bindings_for_container(files_to_bind: Namedlist, base_container_path: str, access: str = "ro") -> tuple[str, dict[str, str]]:
    fns = {}
    for key, path in files_to_bind.items():
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

(input_bindings, input_fns) = bindings_for_container(snakemake.input, base_container_path="/inputs")
(output_bindings, output_fns) = bindings_for_container(snakemake.output, base_container_path="/outputs", access="rw")

alleles = ",".join(args["alleles"])

if peptides := getattr(snakemake.input, "peptides", ""):
    peptides = f"--peptide-fasta {peptides}"
if genes := getattr(snakemake.input, "genes", ""):
    genes = f"--genes-of-interest-file {genes}"

n_threads = min(args["n_threads"], snakemake.threads)

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

export TMPDIR=$(realpath $TMPDIR)

# Re-home pVACtools to avoid conflicts with user's setting (.bashrc, ...)
# Create home on TMPDIR? But then the calling script is gone...
home=$(dirname {snakemake.output})/../home
home=$(realpath $home)
rm -rf $home
mkdir -p $home

rm -rf $(dirname {snakemake.output})
mkdir -p $(dirname {snakemake.output})

cat << __EOF > $home/run_pVACseq.sh
pvacseq run --n-threads {n_threads} \\
    --normal-sample-name {args[normal_sample]} \\
    --iedb-install-directory /opt/iedb \\
    {args[extra_args]} \\
    {peptides} {genes} \\
    {input_fns[vcf]} \\
    {args[tumor_sample]} {alleles} {args[algorithms]} \\
    $(dirname {output_fns[done]})
__EOF
chmod +x $home/run_pVACseq.sh

apptainer exec \
    --home $home --bind $TMPDIR:$TMPDIR:rw \
    {input_bindings} {output_bindings} {snakemake.input[container]} bash $home/run_pVACseq.sh
"""
)

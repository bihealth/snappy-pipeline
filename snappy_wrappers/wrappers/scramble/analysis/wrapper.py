"""CUBI+Snakemake wrapper code for scramble (analysis): Snakemake wrapper.py"""
import os

from snakemake import shell

shell.executable("/bin/bash")

# Define input full path
input_full_path = os.path.join(os.getcwd(), str(snakemake.input))

# Define prefix based on input
prefix = input_full_path.replace("_cluster.txt", "")

shell(
    r"""
set -x

# Pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# Create out dir
mkdir -p $(dirname {snakemake.output.txt})

# Call tool
Rscript --vanilla {snakemake.params.args[rscript]}  --out-name {prefix} \
    --cluster-file {input_full_path} \
    --install-dir {snakemake.config[step_config][mobile_element_insertion][scramble_install_dir]} \
    --mei-refs {snakemake.params.args[mei_refs]} \
    --nCluster {snakemake.params.args[n_cluster]} \
    --mei-score {snakemake.params.args[mei_score]} \
    --indel-score {snakemake.params.args[indel_score]} \
    --poly-a-frac {snakemake.params.args[mei_polya_frac]} \
    --eval-meis
"""
)

# Compute MD5 sums of log and MEI txt.
shell(
    r"""
md5sum {snakemake.log} > {snakemake.log}.md5
md5sum {snakemake.output.txt} > {snakemake.output.txt_md5}
"""
)
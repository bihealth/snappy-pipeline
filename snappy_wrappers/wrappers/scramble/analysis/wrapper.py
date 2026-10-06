"""CUBI+Snakemake wrapper code for scramble (analysis): Snakemake wrapper.py"""

import os

from snakemake import shell

shell.executable("/bin/bash")

args = getattr(snakemake.params, "args", {})

# Define input full path
input_full_path = os.path.join(os.getcwd(), str(snakemake.input))

# Define prefix based on input
prefix = input_full_path.replace("_cluster.txt", "")

# Include user provided MEI Ref if any
mei_ref_argument = ""
if args["mei_refs"]:
    mei_ref_argument = "--mei-refs " + str(args["mei_refs"])

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
scramble.sh  {mei_ref_argument} \
  --ref {args[reference_genome]} \
  --out-name {prefix} \
  --cluster-file {input_full_path} \
  --nCluster {args[n_cluster]} \
  --mei-score {args[mei_score]} \
  --indel-score {args[indel_score]} \
  --poly-a-frac {args[mei_polya_frac]} \
  --eval-meis

# Post-process VCF
bgzip --stdout {snakemake.output.vcf} > {snakemake.output.vcf_gz}
tabix {snakemake.output.vcf_gz}
"""
)

# Compute MD5 sums of log and MEI output
shell(
    r"""
md5sum {snakemake.log} > {snakemake.log}.md5
md5sum {snakemake.output.txt} > {snakemake.output.txt_md5}
md5sum {snakemake.output.vcf_gz} > {snakemake.output.vcf_gz_md5}
md5sum {snakemake.output.vcf_tbi} > {snakemake.output.vcf_tbi_md5}
"""
)

"""CUBI+Snakemake wrapper code for scramble (analysis): Snakemake wrapper.py"""

import os
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

# Define input full path
input_full_path = os.path.join(os.getcwd(), str(snakemake.input.cluster))

# Define prefix based on input
prefix = input_full_path.replace(".cluster.txt", "")

# Include user provided MEI Ref if any
mei_ref_argument = ""
if mei_refs := snakemake.input.get("mei_refs", ""):
    mei_ref_argument = f"--mei-refs {mei_refs}"

ShellWrapper(snakemake).run(
    r"""
# Create out dir
mkdir -p $(dirname {snakemake.output.txt})

# Call tool
scramble.sh  {mei_ref_argument} \
  --ref {snakemake.input.reference_genome} \
  --out-name {prefix} \
  --cluster-file {input_full_path} \
  --nCluster {args[n_cluster]} \
  --mei-score {args[mei_score]} \
  --indel-score {args[indel_score]} \
  --poly-a-frac {args[mei_polya_frac]} \
  --eval-meis

# Scramble names its table {{prefix}}_MEIs.txt, the declared name is {{prefix}}.MEIs.txt
mv {prefix}_MEIs.txt {snakemake.output.txt}

# Post-process VCF
bgzip --stdout {snakemake.output.vcf} > {snakemake.output.vcf_gz}
tabix {snakemake.output.vcf_gz}
"""
)

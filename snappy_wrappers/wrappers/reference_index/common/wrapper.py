from snakemake.shell import shell

shell(
    r"""
set -euo pipefail
set -x

mkdir -p $(dirname {snakemake.log.log})
exec 2> >(tee -a {snakemake.log.log} >&2)

conda info > {snakemake.log.conda_info}
conda list > {snakemake.log.conda_list}

mkdir -p work/reference_index/out
samtools faidx {snakemake.input.reference}
cp -f {snakemake.input.reference}.fai {snakemake.output.reference_fai}
samtools dict -o {snakemake.output.reference_dict} {snakemake.input.reference}
cut -f1,2 {snakemake.output.reference_fai} > {snakemake.output.reference_genome}

md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}

for dst in {snakemake.output.output_links}; do
  src=work/${dst#output/}
  mkdir -p "$(dirname "$dst")"
  ln -snrf "$src" "$dst"
done
"""
)


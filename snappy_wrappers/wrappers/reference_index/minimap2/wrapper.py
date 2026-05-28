from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

extra_args = snakemake.params.args.get("extra_args", "")

shell(
    r"""
set -euo pipefail
set -x

mkdir -p $(dirname {snakemake.log.log})
exec 2> >(tee -a {snakemake.log.log} >&2)

conda info > {snakemake.log.conda_info}
conda list > {snakemake.log.conda_list}

minimap2 -t {snakemake.threads} -d {snakemake.output.mmi} {extra_args} {snakemake.input.reference}

for dst in {snakemake.output.output_links}; do
  src=${{dst/\/output\//\/work\/}}
  mkdir -p "$(dirname "$dst")"
  ln -snrf "$src" "$dst"
done

md5sum {snakemake.log.log} > {snakemake.log.log_md5}
md5sum {snakemake.log.conda_info} > {snakemake.log.conda_info_md5}
md5sum {snakemake.log.conda_list} > {snakemake.log.conda_list_md5}
"""
)


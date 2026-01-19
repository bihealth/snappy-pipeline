from typing import TYPE_CHECKING

from snakemake.shell import shell

if TYPE_CHECKING:
    from snakemake.script import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
min_dp = args["min_dp"]
reference_index_path = args["reference_index_path"]

shell(
    r"""
set -x

# Write files for reproducibility -----------------------------------------------------------------

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

# Also pipe stderr to log file --------------------------------------------------------------------

if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Create auto-cleaned temporary directory
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Run actual tools --------------------------------------------------------------------------------

bcftools query \
    -s {snakemake.wildcards.donor_library_name} \
    -f '%CHROM\t%POS[\t%DP\t%AD]\n' \
    {snakemake.input.vcf} \
| awk -F $'\t' 'BEGIN {{ OFS=FS; prev=0; }}
        {{ if (prev != $1) {{
            printf("variableStep chrom=%s\n", $1, span);
        }} else {{
            if (old2 < $2) {{
                dp = old3;
                split(old4, a, ",");
                rd = a[1];
                if (dp >= {min_dp}) {{
                    printf("%s\t%f\n", old2, (dp - rd) / dp);
                }}
            }}
        }}
        old2=$2;
        old3=$3;
        old4=$4;
        prev=$1;
    }}' \
> $TMPDIR/tmp.wig

cut -f 1-2 {reference_index_path} \
> $TMPDIR/chrom.sizes

wigToBigWig $TMPDIR/tmp.wig $TMPDIR/chrom.sizes {snakemake.output.bw}
pushd $(dirname {snakemake.output.bw})
md5sum $(basename {snakemake.output.bw}) >$(basename {snakemake.output.bw_md5})
popd

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
sleep 1s  # try to wait for log file flush
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
set -x

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

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

WINDOW={snakemake.config[step_config][ngs_mapping][bam_collect_doc][window_length]}

# Compute coverage vcf.gz

maelstrom-core \
    bam-collect-doc \
    --in {snakemake.input.bam} \
    --out {snakemake.output.vcf} \
    --reference {snakemake.config[static_data_config][reference][path]} \
    --window-length $WINDOW

find $(dirname $(dirname {snakemake.output.vcf}))

pushd $(dirname {snakemake.output.vcf})
tabix -f $(basename {snakemake.output.vcf})

md5sum $(basename {snakemake.output.vcf}) >$(basename {snakemake.output.vcf_md5})
md5sum $(basename {snakemake.output.tbi}) >$(basename {snakemake.output.tbi_md5})

# Convert to bigWig file

bcftools query -f '%CHROM\t%POS[\t%CV]\n' $(basename {snakemake.output.vcf}) \
| awk -v span=$WINDOW -F $'\t' 'BEGIN {{ OFS=FS; prev=0; }}
        {{ if (prev != $1) {{
            printf("variableStep chrom=%s span=%d\n", $1, span);
        }} else {{
            printf("%s\t%f\n", old2, old3);
        }}
        old2=$2;
        old3=$3;
        prev=$1;
    }}' \
> $TMPDIR/out.wig
cut -f 1-2 {snakemake.config[static_data_config][reference][path]}.fai \
> $TMPDIR/chrom.sizes

wigToBigWig $TMPDIR/out.wig $TMPDIR/chrom.sizes $(basename {snakemake.output.bw})
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

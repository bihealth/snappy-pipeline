# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log.log}")
set -x
# -----------------------------------------------------------------------------

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

bcftools query \
    -s {snakemake.wildcards.donor_ngs_library} \
    -f '%CHROM\t%POS[\t%DP\t%AD]\n' \
    {snakemake.input.vcf} \
| awk -F $'\t' 'BEGIN {{ OFS=FS; prev=0; }}
        {{ if (prev != $1) {{
            printf("variableStep chrom=%s\n", $1, span);
        }} else {{
            dp = old3;
            split(old4, a, ",");
            rd = a[1];
            if (dp >= {snakemake.config[step_config][variant_calling][baf_file_generation][min_dp]}) {{
                printf("%s\t%f\n", old2, (dp - rd) / dp);
            }}
        }}
        old2=$2;
        old3=$3;
        old4=$4;
        prev=$1;
    }}' \
> $TMPDIR/tmp.wig

cut -f 1-2 {snakemake.config[static_data_config][reference][path]}.fai \
> $TMPDIR/chrom.sizes

wigToBigWig $TMPDIR/tmp.wig $TMPDIR/chrom.sizes {snakemake.output.bw}
pushd $(dirname {snakemake.output.bw})
md5sum $(basename {snakemake.output.bw}) >$(basename {snakemake.output.bw_md5})
"""
)

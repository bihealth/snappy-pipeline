# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
# -----------------------------------------------------------------------------
# Redirect stderr to log file by default and enable printing executed commands
exec 2> >(tee -a "{snakemake.log}")
set -x
# -----------------------------------------------------------------------------

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
            if (dp > 0) {{
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

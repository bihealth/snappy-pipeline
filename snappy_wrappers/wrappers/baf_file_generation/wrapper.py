from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
min_dp = args["min_dp"]
reference_index_path = args["reference_index_path"]

ShellWrapper(snakemake).run(
    r"""
set -x
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
"""
)

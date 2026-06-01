# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

ShellWrapper(snakemake).run(
    r"""
WINDOW={args[window_length]}

# Compute coverage vcf.gz

maelstrom-core \
    bam-collect-doc \
    --in {snakemake.input.bam} \
    --out {snakemake.output.vcf} \
    --reference {snakemake.input.reference} \
    --window-length $WINDOW

find $(dirname $(dirname {snakemake.output.vcf}))

pushd $(dirname {snakemake.output.vcf})
tabix -f $(basename {snakemake.output.vcf})

# Convert coverage to bigWig file

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
> $TMPDIR/out_cov.wig
cut -f 1-2 {snakemake.input.reference}.fai \
> $TMPDIR/chrom.sizes

wigToBigWig $TMPDIR/out_cov.wig $TMPDIR/chrom.sizes $(basename {snakemake.output.cov_bw})

# Convert mapping quality to bigWig file

bcftools query -f '%CHROM\t%POS[\t%MQ]\n' $(basename {snakemake.output.vcf}) \
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
> $TMPDIR/out_mq.wig
cut -f 1-2 {snakemake.input.reference}.fai \
> $TMPDIR/chrom.sizes

wigToBigWig $TMPDIR/out_mq.wig $TMPDIR/chrom.sizes $(basename {snakemake.output.mq_bw})

popd
"""
)


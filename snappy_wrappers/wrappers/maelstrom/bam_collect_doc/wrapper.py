# -*- coding: utf-8 -*-

from snakemake.shell import shell

shell(
    r"""
set -x

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
"""
)

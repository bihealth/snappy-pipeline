from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

ShellWrapper(snakemake).run(
    r"""
set -x

delly merge \
    --outfile $TMPDIR/tmp.bcf \
    {snakemake.input.bcf}

# Some yak-shaving for removing SVs starting at 0 which BCF does not like.
bcftools view \
    -O v \
    $TMPDIR/tmp.bcf \
| awk -F $'\t' '
    BEGIN {{ OFS=FS; }}
    (/^#/ || ($2 != 0)) {{ print $0; }}
    ($2 == 0) {{ $2 = 1; print $0; }}' \
| bcftools view \
    -O b \
    /dev/stdin \
> {snakemake.output.bcf}

tabix -f {snakemake.output.bcf}
"""
)

# -*- coding: utf-8 -*-
"""Wrapper for running Delly2's call merging step
"""

import tempfile

from snakemake.shell import shell

__author__ = "Manuel Holtgrewe"
__email__ = "manuel.holtgrewe@bihealth.de"

# Bad regions on hg19 with too high BND call count.
bad_hg19 = (
    ("2", 33091770, 33092140, "troublesome poly-T"),
    ("2", 33141280, 33141700, "troublesome poly-G"),
    ("12", 66451360, 66451390, "troublesome poly-T"),
    ("12", 66451440, 66451470, "troublesome ALU"),
    ("15", 55218220, 55218290, "troublesome L1"),
    ("15", 77910860, 77910960, "troublesome poly-T+ALU"),
)
reference = snakemake.config["static_data_config"]["reference"]["path"]
if "hg19" in reference or "GRCh37" in reference:
    fix_bnds = r"""
        | awk -F $'\t' '
            /^#/ {{ print $0 }}
            /^[^#]/ {{
                if (($1 == "2" && $2 >= 33091770 && $2 <= 33092140) || ($1 == "2" && $2 >= 33141280 && $2 <= 33141700) || ($1 == "12" && $2 >= 66451360 && $2 <= 66451390) || ($1 == "12" && $2 >= 66451440 && $2 <= 66451470) || ($1 == "15" && $2 >= 55218220 && $2 <= 55218290) || ($1 == "15" && $2 >= 77910860 && $2 <= 77910960)) {{
                    next;
                }} else {{
                    print $0;
                }}
            }}'
    """.strip()
else:
    fix_bnds = ""

with tempfile.NamedTemporaryFile("wt") as tmpf:
    # Write paths to input files into temporary file.
    #
    # cf. https://bitbucket.org/snakemake/snakemake/issues/878
    print("\n".join(snakemake.input), file=tmpf)
    tmpf.flush()
    # Actually run the script.
    shell(
        r"""
    # -----------------------------------------------------------------------------
    # Redirect stderr to log file by default and enable printing executed commands
    exec &> >(tee -a "{snakemake.log}")
    set -x
    # -----------------------------------------------------------------------------

    export LC_ALL=C
    export TMPDIR=$(mktemp -d)
    trap "rm -rf $TMPDIR" EXIT

    mkdir $TMPDIR/cwd

    i=0
    for x in $(cat {tmpf.name}); do
        let "i=$i+1"
        ln -s $(readlink -f $x) $TMPDIR/cwd/$i.bcf
        ln -s $(readlink -f $x).csi $TMPDIR/cwd/$i.bcf.csi
    done

    pushd $TMPDIR/cwd
    delly merge \
        --outfile $TMPDIR/tmp.bcf \
        *.bcf
    popd

    # Some yak-shaving for removing SVs starting at 0 which BCF does not like.
    bcftools view \
        -O v \
        $TMPDIR/tmp.bcf \
    {fix_bnds} \
    | awk -F $'\t' '
        BEGIN {{ OFS=FS; }}
        (/^#/ || ($2 != 0)) {{ print $0; }}
        ($2 == 0) {{ $2 = 1; print $0; }}' \
    | bcftools view \
        -O b \
        /dev/stdin \
    > {snakemake.output.bcf}

    tabix -f {snakemake.output.bcf}

    pushd $(dirname {snakemake.output.bcf})
    md5sum $(basename {snakemake.output.bcf}) >$(basename {snakemake.output.bcf}).md5
    md5sum $(basename {snakemake.output.bcf}).csi >$(basename {snakemake.output.bcf}).csi.md5
    """
    )

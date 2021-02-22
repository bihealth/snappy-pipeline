# -*- coding: utf-8 -*-
"""Wrapper for running CopywriteR"""

from snakemake import shell

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bihealth.de>"

shell.executable("/bin/bash")

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT

# Also pipe stderr to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec 2> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# -------------------------------------------------------------------------------------------------
# Write helper script and call R
#
cat <<"EOF" > $TMPDIR/run_copywriter.R
library( CopywriteR )
preCopywriteR( output.folder='.',
    bin.size={snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][bin_size]},
    ref.genome='{snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][genome]}',
    prefix='{snakemake.config[step_config][somatic_targeted_seq_cnv_calling][copywriter][prefix]}')

warnings()
sessionInfo()
EOF

pushd $TMPDIR
Rscript --vanilla run_copywriter.R
popd

# -------------------------------------------------------------------------------------------------
# Move out output
#
# copywriter puts the files into a directory named $genome_release_$binsize/1000/kb, e.g. hg19_20kb

mv $TMPDIR/*/GC_mappability.rda {snakemake.output.gc}
mv $TMPDIR/*/blacklist.rda {snakemake.output.blacklist}

# md5 sums
out_dir=$(readlink -f $(dirname {snakemake.output.gc}))
pushd $out_dir && \
    for i in *.rda; do
        md5sum $i > ${{i}}.md5
    done
    popd
"""
)

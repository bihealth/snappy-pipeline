# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing PureCN panel of normals"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

genomicsDB = getattr(snakemake.input, "genomicsdb", "")
if genomicsDB:
    # The archive is extracted in $TMPDIR, so its path must not be relative
    genomicsDB = os.path.abspath(genomicsDB)

ShellWrapper(snakemake).run(
    r"""
outdir=$TMPDIR/out
mkdir -p $outdir

mkdir $TMPDIR/extra
echo "{snakemake.input.normals}" | tr " " "\n" > $TMPDIR/extra/filenames.txt

# Extract Mutect2 genomicsDB when present
normal_panel=""
if [[ -n "{genomicsDB}" ]]
then
    pushd $TMPDIR ; tar -zxvf {genomicsDB} ; popd
    normal_panel=" --normal-panel $TMPDIR/pon_db "
fi

# Create panel
Rscript /opt/PureCN/NormalDB.R \
    --out-dir $outdir \
    --coverage-files $TMPDIR/extra/filenames.txt $normal_panel \
    --genome {config[genome_name]} --assay {config[enrichment_kit_name]} \
    {extra_arguments}

# Move output to destination
mv $outdir/normalDB_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.db}
mv $outdir/mapping_bias_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.mapbias}
mv $outdir/mapping_bias_hq_sites_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.hq}
mv $outdir/low_coverage_targets_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.lowcov}
mv $outdir/interval_weights_{config[enrichment_kit_name]}_{config[genome_name]}.png {snakemake.output.plot}
"""
)

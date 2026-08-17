# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing PureCN panel of normals"""

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

genomicsDB = getattr(snakemake.input, "genomicsdb", "")

ShellWrapper(snakemake).run(
    r"""
outdir=$TMPDIR/out
mkdir -p $outdir

mkdir $TMPDIR/extra
echo "{snakemake.input.normals}" | tr " " "\n" > $TMPDIR/extra/filenames.txt

# Extract Mutect2 genomicsDB when present
if [[ -n "{genomicsDB}" ]]
then
    pushd $TMPDIR ; tar -zxvf {genomicsDB} ; popd
    normal_panel=" --normal-panel /pon_db "
else
    mkdir $TMPDIR/pon_db
    normal_panel=""
fi

# Create panel
cmd="/usr/local/bin/Rscript /opt/PureCN/NormalDB.R \
    --out-dir /output \
    --coverage-files /extra/filenames.txt $normal_panel \
    --genome {config[genome_name]} --assay {config[enrichment_kit_name]}
"
apptainer exec --home $PWD -B $outdir:/output -B $TMPDIR/pon_db:/pon_db:ro -B $TMPDIR/extra:/extra:ro {snakemake.input.container} $cmd

# Move output to destination
mv $outdir/normalDB_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.db}
mv $outdir/mapping_bias_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.mapbias}
mv $outdir/mapping_bias_hq_sites_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.hq}
mv $outdir/low_coverage_targets_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.lowcov}
mv $outdir/interval_weights_{config[enrichment_kit_name]}_{config[genome_name]}.png {snakemake.output.plot}
"""
)

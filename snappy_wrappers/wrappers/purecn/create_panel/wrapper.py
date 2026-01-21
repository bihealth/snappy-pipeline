# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for preparing PureCN panel of normals"""

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

if "path_genomicsDB" in config.keys() and config["path_genomicsDB"]:
    genomicsDB = config["path_genomicsDB"]
else:
    genomicsDB = ""

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

# MD5 checksums without dirname
md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo $checksum
}}

outdir=$tmpdir/out
mkdir -p $outdir

mkdir $tmpdir/extra
echo "{snakemake.input.normals}" | tr " " "\n" > $tmpdir/extra/filenames.txt

# Extract Mutect2 genomicsDB when present
if [[ -n "{genomicsDB}" ]]
then
    pushd $tmpdir ; tar -zxvf {genomicsDB} ; popd
    normal_panel=" --normal-panel /pon_db "
else
    mkdir $tmpdir/pon_db
    normal_panel=""
fi

# Create panel
cmd="/usr/local/bin/Rscript /opt/PureCN/NormalDB.R \
    --out-dir /output \
    --coverage-files /extra/filenames.txt $normal_panel \
    --genome {config[genome_name]} --assay {config[enrichment_kit_name]}
"
apptainer exec --home $PWD -B $outdir:/output -B $tmpdir/pon_db:/pon_db:ro -B $tmpdir/extra:/extra:ro {snakemake.input.container} $cmd

# Move output to destination
mv $outdir/normalDB_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.db}
mv $outdir/mapping_bias_{config[enrichment_kit_name]}_{config[genome_name]}.rds {snakemake.output.mapbias}
mv $outdir/mapping_bias_hq_sites_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.hq}
mv $outdir/low_coverage_targets_{config[enrichment_kit_name]}_{config[genome_name]}.bed {snakemake.output.lowcov}
mv $outdir/interval_weights_{config[enrichment_kit_name]}_{config[genome_name]}.png {snakemake.output.plot}

# MD5 checksum for main result only
md5 {snakemake.output.db} > {snakemake.output.db}.md5
md5 {snakemake.output.mapbias} > {snakemake.output.mapbias}.md5
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

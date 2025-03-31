# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for manta/strelka2: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

this_file = __file__

args = getattr(snakemake.params, "args", {})

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        exec &> >(tee -a "{snakemake.log}" >&2)
    else
        rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
        echo "No tty, logging disabled" >"{snakemake.log}"
    fi
fi

# TODO: add through shell.prefix
export TMPDIR=/fast/users/$USER/scratch/tmp

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/manta
mkdir -p $TMPDIR/strelka
mkdir -p $TMPDIR/out
mkdir -p $TMPDIR/log

# Extract the base filename (typically <mapper>.strelka.<tumor_library>)
out_base=$(basename "{snakemake.output.vcf}" .vcf.gz)

# Configure & run manta to get indels candidates
configManta.py \
    --normalBam "{snakemake.input.normal_bam}" \
    --tumorBam  "{snakemake.input.tumor_bam}"  \
    --referenceFasta "{args[reference]}" \
    --runDir $TMPDIR/manta

$TMPDIR/manta/runWorkflow.py -m local -j 8

# Use target bed file if present
cmd=""
if [[ "X{args[path_target_regions]}" != "X" ]]
then
    if [[ -r "{args[path_target_regions]}" ]]
    then
        cmd=' --exome --callRegions "{args[path_target_regions]}" '
    fi
fi

# Configure & run strelka2
configureStrelkaSomaticWorkflow.py \
    --normalBam "{snakemake.input.normal_bam}" \
    --tumorBam  "{snakemake.input.tumor_bam}"  \
    --referenceFasta "{args[reference]}" \
    --indelCandidates "$TMPDIR/manta/results/variants/candidateSmallIndels.vcf.gz" \
    $cmd --outputCallableRegions \
    --runDir $TMPDIR/strelka

$TMPDIR/strelka/runWorkflow.py -m local -j 8

# Post-processing
# 1. Merge SNVs & indels
# 2. Add library names for normal & tumor column in output vcfs
# 3. Filter only PASS variants

cat <<"EOF" >$TMPDIR/samples.txt
{snakemake.params.tumor_lib_name}
{snakemake.params.normal_lib_name}
EOF

bcftools concat \
    --allow-overlaps --rm-dups all \
    $TMPDIR/strelka/results/variants/somatic.snvs.vcf.gz \
    $TMPDIR/strelka/results/variants/somatic.indels.vcf.gz \
    | bcftools reheader \
    --samples $TMPDIR/samples.txt \
    --output $TMPDIR/out/${{out_base}}.full.vcf
bgzip $TMPDIR/out/${{out_base}}.full.vcf
tabix $TMPDIR/out/${{out_base}}.full.vcf.gz

bcftools filter \
    --include 'FILTER=="PASS"' \
    --output $TMPDIR/out/${{out_base}}.vcf.gz --output-type z \
    $TMPDIR/out/${{out_base}}.full.vcf.gz
tabix $TMPDIR/out/${{out_base}}.vcf.gz

# Move extra files
mv $TMPDIR/strelka/results/regions/somatic.callable.regions.bed.gz $TMPDIR/out/${{out_base}}.bed.gz
mv $TMPDIR/strelka/results/regions/somatic.callable.regions.bed.gz.tbi $TMPDIR/out/${{out_base}}.bed.gz.tbi
mv $TMPDIR/strelka/results/stats/runStats.tsv $TMPDIR/out/${{out_base}}.tsv
mv $TMPDIR/strelka/results/stats/runStats.xml $TMPDIR/out/${{out_base}}.xml

# MD5 checksums & mode files to final destination
pushd $TMPDIR/out && \
    for f in ${{out_base}}.*; do \
        md5sum $f >$f.md5; \
    done && \
    popd

mv $TMPDIR/out/${{out_base}}.* $(dirname {snakemake.output.vcf})

# Additional logging for transparency & reproducibility
# Logging: Save a copy this wrapper (with the pickle details in the header)
cp {this_file} $(dirname {snakemake.log.log})/wrapper.py

# Logging: Save a permanent copy of the environment file used
cp $(dirname {this_file})/environment.yaml $(dirname {snakemake.log.log})/environment_wrapper.yaml

pushd $(dirname {snakemake.log.log})
for f in $(ls)
do
    md5sum $f > $f.md5
done
popd
"""
)

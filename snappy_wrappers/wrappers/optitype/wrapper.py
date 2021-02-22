# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Optitype: Snakemake wrapper.py
"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = snakemake.params.args["input"]["reads_left"]
reads_right = snakemake.params.args["input"].get("reads_right", "")

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

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

if [[ -z "{reads_right}" ]]; then
    paired=0
else
    paired=1
fi

# First alignment step (filter to candidate HLA reads) --------------------------------------------

data_dir=$(dirname $(readlink -f $(which OptiTypePipeline.py)))/data
seq_type={snakemake.params.args[seq_type]}
refname=hla_reference_$seq_type.fasta
reference=$data_dir/$refname

mkdir -p reference $TMPDIR/yara_index
cp $reference $TMPDIR/yara_index
yara_indexer -o $TMPDIR/yara_index/$refname $TMPDIR/yara_index/$refname

yara_mapper \
    -t {snakemake.config[step_config][hla_typing][optitype][num_mapping_threads]} \
    --error-rate 5 \
    --output-format sam \
    $TMPDIR/yara_index/$refname \
    <(zcat --force {reads_left}) \
| samtools fastq -F 4 - \
> $TMPDIR/tmp.d/reads_left.fastq

if [[ $paired -eq 1 ]]; then
    yara_mapper \
        -t {snakemake.config[step_config][hla_typing][optitype][num_mapping_threads]} \
        --error-rate 5 \
        --output-format sam \
        $TMPDIR/yara_index/$refname \
        <(zcat --force {reads_right}) \
    | samtools fastq -F 4 - \
    > $TMPDIR/tmp.d/reads_right.fastq
fi

cat $TMPDIR/tmp.d/reads_*.fastq | \
seqtk sample - {snakemake.config[step_config][hla_typing][optitype][max_reads]} \
> $TMPDIR/tmp.d/reads_sampled.fastq

wc -l $TMPDIR/tmp.d/reads_*.fastq

# Run OptiType ------------------------------------------------------------------------------------

cat <<"EOF" >$TMPDIR/tmp.d/optitype.ini
[mapping]
razers3=razers3
threads={snakemake.config[step_config][hla_typing][optitype][num_mapping_threads]}

[ilp]
solver=glpk
threads=1

[behavior]
deletebam=true
unpaired_weight=0
use_discordant=false
EOF

mkdir -p $TMPDIR/out.tmp
OptiTypePipeline.py \
    --config $TMPDIR/tmp.d/optitype.ini \
    --verbose \
    --input $TMPDIR/tmp.d/reads_sampled.fastq \
    --$seq_type \
    --outdir $TMPDIR/out.tmp

pushd $TMPDIR/out.tmp/*
for f in *; do
    test -f $f && md5sum $f >$f.md5
done
popd

# move files to the output directory
prefix=$(basename $(ls $TMPDIR/out.tmp | head -n 1))
mv $TMPDIR/out.tmp/${{prefix}}/${{prefix}}_coverage_plot.pdf {snakemake.output.cov_pdf}
mv $TMPDIR/out.tmp/${{prefix}}/${{prefix}}_coverage_plot.pdf.md5 {snakemake.output.cov_pdf}.md5
mv $TMPDIR/out.tmp/${{prefix}}/${{prefix}}_result.tsv {snakemake.output.tsv}
mv $TMPDIR/out.tmp/${{prefix}}/${{prefix}}_result.tsv.md5 {snakemake.output.tsv}.md5

# create final .txt file
tail -n +2 {snakemake.output.tsv} \
    | cut -f 2-7 \
    | tr '\t' '\n' \
    | sort \
    > {snakemake.output.txt}
md5sum {snakemake.output.txt} > {snakemake.output.txt_md5}
"""
)

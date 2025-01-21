# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for Optitype: Snakemake wrapper.py"""

from snakemake import shell

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
reads_left = args["input"]["reads_left"]
reads_right = args["input"].get("reads_right", "")

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Setup auto-cleaned TMPDIR
export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" EXIT
mkdir -p $TMPDIR/tmp.d

# Also pipe stderr to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec 2> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_list}.md5
md5sum {snakemake.log.conda_info} | sed -re "s/  (\.?.+\/)([^\/]+)$/  \2/" > {snakemake.log.conda_info}.md5

if [[ -z "{reads_right}" ]]; then
    paired=0
else
    paired=1
fi

# First alignment step (filter to candidate HLA reads) --------------------------------------------

data_dir=$(dirname $(readlink -f $(which OptiTypePipeline.py)))/data
seq_type={args[seq_type]}
refname=hla_reference_$seq_type.fasta
reference=$data_dir/$refname

mkdir -p reference $TMPDIR/yara_index
cp $reference $TMPDIR/yara_index
yara_indexer -o $TMPDIR/yara_index/$refname $TMPDIR/yara_index/$refname

yara_mapper \
    -t {args[num_mapping_threads]} \
    --error-rate 5 \
    --output-format sam \
    $TMPDIR/yara_index/$refname \
    <(zcat --force {reads_left}) \
| samtools fastq -F 4 - \
> $TMPDIR/tmp.d/reads_left.fastq

if [[ $paired -eq 1 ]]; then
    yara_mapper \
        -t {args[num_mapping_threads]} \
        --error-rate 5 \
        --output-format sam \
        $TMPDIR/yara_index/$refname \
        <(zcat --force {reads_right}) \
    | samtools fastq -F 4 - \
    > $TMPDIR/tmp.d/reads_right.fastq
fi

cat $TMPDIR/tmp.d/reads_*.fastq | \
seqtk sample - {args[max_reads]} \
> $TMPDIR/tmp.d/reads_sampled.fastq

wc -l $TMPDIR/tmp.d/reads_*.fastq

# Run OptiType ------------------------------------------------------------------------------------

cat <<"EOF" >$TMPDIR/tmp.d/optitype.ini
[mapping]
razers3=razers3
threads={args[num_mapping_threads]}

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

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

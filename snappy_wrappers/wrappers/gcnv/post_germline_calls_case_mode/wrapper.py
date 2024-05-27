# -*- coding: utf-8 -*-

from snakemake.shell import shell

paths_calls = " ".join(snakemake.input.calls)
paths_models = " ".join(snakemake.params.args["model"])

shell(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"

itv_vcf={snakemake.output.itv_vcf}
seg_vcf={snakemake.output.seg_vcf}

sample_index=-1
for path in $(dirname {snakemake.input.calls[0]})/cnv_calls-calls/SAMPLE_*; do
    if [[ "$(cat $path/sample_name.txt)" == "{snakemake.wildcards.library_name}" ]]; then
        sample_index=$(basename $path | sed -e 's/SAMPLE_//')
        break
    fi
done
if [[ $sample_index == -1 ]]; then
    >&2 echo "Could not determine sample index."
    exit 1
fi

# Get contig name style; SAMPLE_0  is guaranteed to exist
tail -n +2 $(dirname {snakemake.input.ploidy})/ploidy-calls/SAMPLE_0/contig_ploidy.tsv | grep "^chr[0-9XY]{{1,2}}\s[0-9]\s[0-9.]+$" > /dev/null && true
exit_status=$?
if [[ $exit_status == 0 ]]; then
    STYLE="chr"
else
    STYLE=""
fi

gatk PostprocessGermlineCNVCalls \
    $(for x in {paths_calls}; do echo --calls-shard-path $(dirname $x)/cnv_calls-calls; done) \
    $(for x in {paths_models}; do echo --model-shard-path $x; done) \
    --contig-ploidy-calls $(dirname {snakemake.input.ploidy})/ploidy-calls \
    --sample-index $sample_index \
    --autosomal-ref-copy-number 2 \
    --allosomal-contig ${{STYLE}}X \
    --allosomal-contig ${{STYLE}}Y \
    --output-genotyped-intervals ${{itv_vcf%.gz}} \
    --output-genotyped-segments ${{seg_vcf%.gz}} \
    --output-denoised-copy-ratios {snakemake.output.ratio_tsv}

perl -p -i -e 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' ${{itv_vcf%.gz}}
perl -p -i -e 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' ${{seg_vcf%.gz}}

bgzip ${{itv_vcf%.gz}}
tabix -f $itv_vcf
bgzip ${{seg_vcf%.gz}}
tabix -f $seg_vcf

for x in $itv_vcf $itv_vcf.tbi $seg_vcf $seg_vcf.tbi {snakemake.output.ratio_tsv}; do
    pushd $(dirname $x)
    md5sum $(basename $x) >$(basename $x).md5
    popd
done
"""
)

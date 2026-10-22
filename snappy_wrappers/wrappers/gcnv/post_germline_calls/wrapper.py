# -*- coding: utf-8 -*-

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
allosomal_contigs = " ".join(args["allosomal_contigs"])
# In cohort mode, each shard directory of call_cnvs holds the calls and the model
paths_shards = " ".join(snakemake.input.calls)

ShellWrapper(snakemake).run(
    r"""
set -x

export TMPDIR=$(mktemp -d)
trap "rm -rf $TMPDIR" ERR EXIT

export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

itv_vcf={snakemake.output.itv_vcf}
seg_vcf={snakemake.output.seg_vcf}

sample_index=-1
for path in {snakemake.input.calls[0]}/cnv_calls-calls/SAMPLE_*; do
    if [[ "$(cat $path/sample_name.txt)" == "{snakemake.wildcards.library_name}" ]]; then
        sample_index=$(basename $path | sed -e 's/SAMPLE_//')
        break
    fi
done

# Get contig name style
# Note: the sample_index in '{snakemake.input.ploidy}' and '{snakemake.input.calls[0]}' are NOT guaranteed to be the same
# -> Sample_0 should always exist and the contig names (which we care for) are the same
tail -n +2 {snakemake.input.ploidy}/ploidy-calls/SAMPLE_0/contig_ploidy.tsv | egrep "^chr[0-9XY]{{1,2}}\s[0-9]\s[0-9.]+" > /dev/null && true
exit_status=$?
if [[ $exit_status == 0 ]]; then
    STYLE="chr"
else
    STYLE=""
fi

gatk PostprocessGermlineCNVCalls \
    $(for x in {paths_shards}; do echo --calls-shard-path $x/cnv_calls-calls; done) \
    $(for x in {paths_shards}; do echo --model-shard-path $x/cnv_calls-model; done) \
    --contig-ploidy-calls {snakemake.input.ploidy}/ploidy-calls \
    --sample-index $sample_index \
    --autosomal-ref-copy-number {args[autosomal_ref_copy_number]} \
    $(for contig in {allosomal_contigs}; do echo --allosomal-contig ${{STYLE}}$contig; done) \
    --output-genotyped-intervals ${{itv_vcf%.gz}} \
    --output-genotyped-segments ${{seg_vcf%.gz}} \
    --output-denoised-copy-ratios {snakemake.output.ratio_tsv}

perl -p -i -e 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' ${{itv_vcf%.gz}}
perl -p -i -e 's/ID=GT,Number=1,Type=Integer/ID=GT,Number=1,Type=String/g' ${{seg_vcf%.gz}}

bgzip ${{itv_vcf%.gz}}
tabix -f $itv_vcf
bgzip ${{seg_vcf%.gz}}
tabix -f $seg_vcf
"""
)

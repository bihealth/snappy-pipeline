from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

bams = " --bam ".join(snakemake.input.bam)

ShellWrapper(snakemake).run(
    r"""
workdir=$TMPDIR/manta

# configManta.py bails out if the run directory exists, so use a new one in the job's TMPDIR
configManta.py \
    --referenceFasta {snakemake.input.reference} \
    --runDir $workdir \
    {extra_arguments} \
    --bam {bams}

perl -p -i -e 's/isEmail = .*/isEmail = False/g' $workdir/runWorkflow.py

python2 $workdir/runWorkflow.py \
    --jobs {snakemake.threads}

cp $workdir/results/variants/diploidSV.vcf.gz {snakemake.output.vcf}
cp $workdir/results/variants/diploidSV.vcf.gz.tbi {snakemake.output.vcf_tbi}
"""
)

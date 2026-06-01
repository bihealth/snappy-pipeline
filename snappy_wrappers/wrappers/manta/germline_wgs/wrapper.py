from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

args = getattr(snakemake.params, "args", {})

bams = " --bam ".join(snakemake.input.bam)

ShellWrapper(snakemake).run(
    r"""
basedir=$(dirname $(dirname {snakemake.output.vcf}))
workdir=$basedir/work
outdir=$basedir/out

# Ensure the working directory is removed, configManta.py will bail out if it already exists
trap "rm -rf \"$workdir\"" EXIT
# Clear out $outdir, there may be some old files remaining that are not governed by Snakemake
rm -rf $outdir/* $workdir/*

configManta.py \
    --referenceFasta {snakemake.input.reference} \
    --runDir $workdir \
    --bam {bams}

perl -p -i -e 's/isEmail = .*/isEmail = False/g' $workdir/runWorkflow.py

python2 $workdir/runWorkflow.py \
    --jobs {snakemake.threads}

cp -ra $workdir/results $outdir
rm -rf $workdir

pushd $outdir
ln -sr results/variants/diploidSV.vcf.gz $(basename {snakemake.output.vcf})
ln -sr results/variants/diploidSV.vcf.gz.tbi $(basename {snakemake.output.vcf_tbi})
ln -sr results/variants/candidateSV.vcf.gz \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz
ln -sr results/variants/candidateSV.vcf.gz.tbi \
    $(basename {snakemake.output.vcf} .vcf.gz).candidates.vcf.gz.tbi
popd
"""
)

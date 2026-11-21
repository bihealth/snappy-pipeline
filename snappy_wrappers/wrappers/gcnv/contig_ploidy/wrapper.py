from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = snakemake.params.args
paths_tsv = " ".join(snakemake.input.tsv)

# Ploidy priors as tab separated values, `echo -e` expands the tabs
autosome_prior = r"\t".join(f"{p:g}" for p in args["autosome_ploidy_prior"])
x_prior = r"\t".join(f"{p:g}" for p in args["x_ploidy_prior"])
y_prior = r"\t".join(f"{p:g}" for p in args["y_ploidy_prior"])

# Add interval block list for PAR regions if configured.
par_intervals = snakemake.input.get("par_intervals", "")
if par_intervals:
    par_args = f"-XL {par_intervals}"
else:
    par_args = ""

ShellWrapper(snakemake).run(
    r"""
export THEANO_FLAGS="base_compiledir=$TMPDIR/theano_compile_dir"
export PYTENSOR_FLAGS="base_compiledir=$TMPDIR/pytensor_compile_dir"

# Get contig name style
egrep "^@SQ\s+SN:chr[0-9XY]{{1,2}}\s+" {snakemake.input.interval_list} > /dev/null && true
exit_status=$?
if [[ $exit_status == 0 ]]; then
    STYLE="chr"
else
    STYLE=""
fi

PRIORS=$TMPDIR/ploidy_priors.tsv
echo -e "CONTIG_NAME\tPLOIDY_PRIOR_0\tPLOIDY_PRIOR_1\tPLOIDY_PRIOR_2\tPLOIDY_PRIOR_3" \
    > $PRIORS
for i in {{1..22}}; do
    echo -e "$STYLE$i\t{autosome_prior}" >> $PRIORS
done
echo -e "${{STYLE}}X\t{x_prior}" >> $PRIORS
echo -e "${{STYLE}}Y\t{y_prior}" >> $PRIORS

mkdir -p {snakemake.output.ploidy}

gatk DetermineGermlineContigPloidy \
    -L {snakemake.input.interval_list} \
    {par_args} \
    --interval-merging-rule {args[interval_merging_rule]} \
    $(for tsv in {paths_tsv}; do echo -I $tsv; done) \
    --contig-ploidy-priors $PRIORS \
    --output {snakemake.output.ploidy} \
    --output-prefix ploidy

# The ploidy calls of X and Y follow the PED files
python $(dirname {__file__})/../../../tools/gcnv_ploidy.py \
    {snakemake.output.ploidy} \
    {snakemake.input.ped}
"""
)

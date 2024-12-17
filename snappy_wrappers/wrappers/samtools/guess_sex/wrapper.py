# -*- coding: utf-8 -*-
"""Wrapper for guessing sex from coverage of autosomes, X & Y"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

cmd=r"""
stat=$TMPDIR/idx_stats.txt

samtools idxstats {snakemake.input.bam} > $stat

aut_len=$(grep -E $'^[0-9]+\t' $stat | cut -f 2 | paste -sd+ | bc)
aut_count=$(grep -E $'^[0-9]+\t' $stat | cut -f 3 | paste -sd+ | bc)
aut_rate=$(echo "$aut_count / $aut_len" | bc -l)

X_len=$(awk -F '\t' '$1 ~ /^(chr)?X$/ {print $2}' $stat)
X_count=$(awk -F '\t' '$1 ~ /^(chr)?X$/ {print $3}' $stat)
x_rate=$(echo "$X_count / $X_len" | bc -l)

Y_len=$(awk -F '\t' '$1 ~ /^(chr)?Y$/ {print $2}' $stat)
Y_count=$(awk -F '\t' '$1 ~ /^(chr)?Y$/ {print $3}' $stat)
Y_rate=$(echo "$Y_count / $Y_len" | bc -l)

cat << __EOF | tr ':' '\t' > {snakemake.output.table}
Name:Length:Count:Rate
Autosomes:$aut_len:$aut_count:$aut_rate
chrX:$X_len:$X_count:$X_rate
chrY:$Y_len:$Y_count:$Y_rate
__EOF

X_ratio=$(echo "2 * $X_rate / $aut_rate" | bc -l)
Y_ratio=$(echo "2 * $Y_rate / $aut_rate" | bc -l)

decision="ambiguous"
if [[ $X_ratio -gt {args[min_X_female]} ]] && [[ $Y_ratio -lt {args[max_Y_female]} ]]
then
    decision="female"
fi
if [[ $X_ratio -lt {args[max_X_male]} ]] && [[ $Y_ratio -gt {args[min_Y_male]} ]]
then
    decision="male"
fi

echo $decision > {snakemake.output.decision}
""".format(
    snakemake=snakemake,
    args=args,
)

ShellWrapper(snakemake).run(cmd)
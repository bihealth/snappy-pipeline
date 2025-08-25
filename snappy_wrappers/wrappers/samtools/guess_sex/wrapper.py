# -*- coding: utf-8 -*-
"""Wrapper for guessing sex from coverage of autosomes, X & Y"""

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

import os
import sys

# The following is required for being able to import snappy_wrappers modules
# inside wrappers. When the wrappers have their own python environment, messing
# with the path is necessary.
base_dir = os.path.normpath(os.path.dirname(__file__))
while os.path.basename(base_dir) != "snappy_wrappers":
    base_dir = os.path.dirname(base_dir)
sys.path.insert(0, os.path.dirname(base_dir))

from snappy_wrappers.snappy_wrapper import ShellWrapper

args = getattr(snakemake.params, "args", {})

cmd=r"""
stat=$TMPDIR/idx_stats.txt

samtools idxstats {snakemake.input.bam} > $stat

aut_len=$(grep -E $'^(chr)?[0-9]+\t' $stat | cut -f 2 | paste -sd+ | bc)
aut_count=$(grep -E $'^(chr)?[0-9]+\t' $stat | cut -f 3 | paste -sd+ | bc)
aut_rate=$(echo "$aut_count / $aut_len" | bc -l)

X_len=$(awk -F '\t' '$1 ~ /^(chr)?X$/ {{print $2}}' $stat)
X_count=$(awk -F '\t' '$1 ~ /^(chr)?X$/ {{print $3}}' $stat)
X_rate=$(echo "$X_count / $X_len" | bc -l)

Y_len=$(awk -F '\t' '$1 ~ /^(chr)?Y$/ {{print $2}}' $stat)
Y_count=$(awk -F '\t' '$1 ~ /^(chr)?Y$/ {{print $3}}' $stat)
Y_rate=$(echo "$Y_count / $Y_len" | bc -l)

X_ratio=$(echo "2 * $X_rate / $aut_rate" | bc -l)
Y_ratio=$(echo "2 * $Y_rate / $aut_rate" | bc -l)

cat << __EOF | tr ':' '\t' > {snakemake.output.table}
Name:Length:Count:Rate:Ratio
Autosomes:$aut_len:$aut_count:$aut_rate:2
chrX:$X_len:$X_count:$X_rate:$X_ratio
chrY:$Y_len:$Y_count:$Y_rate:$Y_ratio
__EOF

# Bash doesn't understand floats
gt() {{
    status=$(echo "$1 > $2" | bc -l)
    test $status -eq 1
}}
lt() {{
    status=$(echo "$1 < $2" | bc -l)
    test $status -eq 1
}}

decision="ambiguous"
if $(gt $X_ratio {args[min_X_female]}) && $(lt $Y_ratio {args[max_Y_female]})
then
    decision="female"
fi
if $(lt $X_ratio {args[max_X_male]}) && $(gt $Y_ratio {args[min_Y_male]})
then
    decision="male"
fi

echo $decision > {snakemake.output.decision}
""".format(
    snakemake=snakemake,
    args=args,
)

ShellWrapper(snakemake).run(cmd)

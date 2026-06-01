# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for infer_experiment: Snakemake wrapper.py"""

import hashlib
from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]
with open(config["path_exon_bed"], "rb") as inputf:
    bed_file_md5 = hashlib.md5(inputf.read()).hexdigest()

ShellWrapper(snakemake).run(
    r"""
# ----- Run rseqc to infer strandedness
infer_experiment.py \
    -r "{config[path_exon_bed]}" \
    -i "{snakemake.input.bam}" \
    > "{snakemake.output.tsv}"

# ----- Parse rseqc infer_experiment output
pattern="^This is (Pair|Single)End Data$"
endedness=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

pattern="^Fraction of reads failed to determine: *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
failed=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")

if [[ "$endedness" = "Single" ]]
then
    pattern="^Fraction of reads explained by \"\\+\\+,\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    pattern="^Fraction of reads explained by \"\\+\\-,\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
else
    pattern="^Fraction of reads explained by \"1\\+\\+,1\\-\\-,2\\+\\-,2\\-\\+\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    forward=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
    pattern="^Fraction of reads explained by \"1\\+\\-,1\\-\\+,2\\+\\+,2\\-\\-\": *([+-]?([0-9]+(\\.[0-9]*)?|\\.[0-9]+)([EeDd][+-]?[0-9]+)?)$"
    reverse=$(grep -E "$pattern" {snakemake.output.tsv} | sed -E "s/$pattern/\1/")
fi

forward=$(echo $forward | tr '[Dde]' 'E')
reverse=$(echo $reverse | tr '[Dde]' 'E')

# ----- Infer protocol strandedness
infer=0
if [ $(echo "$forward > {config[threshold]}" | bc -l) -gt 0 ]
then
    infer=1
fi
if [ $(echo "$reverse > {config[threshold]}" | bc -l) -gt 0 ]
then
    infer=2
fi

decision="{config[strand]}"
if [[ $decision -eq -1 ]]
then
    decision=$infer
fi

# ----- Write outputs
cat << __EOF > {snakemake.output.decision}
{{
    "library_name": "{args[library_name]}",
    "bed_path": "{config[path_exon_bed]}",
    "bed_file_md5": "{bed_file_md5}",
    "bam_path": "{snakemake.input.bam}",
    "strand_from_user": "{config[strand]},
    "strand_from_infer": "$decision",
    "decision_threshold": {config[threshold]},
    "endedness": "$endedness",
    "fraction_forward": $forward,
    "fraction_reverse": $reverse,
    "fraction_failed": $failed,
    "decision": "$decision"
}}
__EOF
"""
)


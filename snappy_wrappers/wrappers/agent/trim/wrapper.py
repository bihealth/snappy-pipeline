"""Trimming with AGeNT, which moves the molecular barcodes into the read names, as mbcs ran it"""

from typing import TYPE_CHECKING

from snappy_wrappers.snappy_wrapper import ShellWrapper

if TYPE_CHECKING:
    from snakemake.iocontainers import snakemake

args = getattr(snakemake.params, "args", {})

# The reads are inputs. If another task wrote them, the input is that task's out directory and
# params list the files. Params also name the trimmed files, in the order of the inputs.
reads = args.get("input", {})
input_path_left = list(snakemake.input.get("reads_left") or reads["reads_left"])
input_path_right = list(snakemake.input.get("reads_right") or reads.get("reads_right", []))
output_left = [f"{x['relative_path']}/{x['filename']}" for x in args["output_left"]]
output_right = [f"{x['relative_path']}/{x['filename']}" for x in args["output_right"]]
assert len(input_path_left) == len(output_left)
assert len(input_path_right) == len(output_right)
config = args["config"]

ShellWrapper(snakemake).run(
    r"""
declare -a reads_left=({input_path_left})
declare -a reads_right=({input_path_right})
declare -a output_left=({output_left})
declare -a output_right=({output_right})

# AGeNT writes <prefix>_R1.fastq.gz and <prefix>_R2.fastq.gz; the trimmed files keep the names
# of the input files
mkdir -p $TMPDIR/trimmed
for ((i = 0; i < ${{#reads_left[@]}}; i++)); do
    fq2=
    if [[ ${{#reads_right[@]}} -gt 0 ]]; then
        fq2="-fq2 ${{reads_right[$i]}}"
    fi
    java {config[java_options]} -Djava.io.tmpdir=$TMPDIR -jar {snakemake.input.path} \
        -{config[lib_prep_type]} {extra_arguments} \
        -out $TMPDIR/trimmed/$i \
        -fq1 ${{reads_left[$i]}} $fq2

    mkdir -p $(dirname {snakemake.output.out}/${{output_left[$i]}})
    mv $TMPDIR/trimmed/${{i}}_R1.fastq.gz {snakemake.output.out}/${{output_left[$i]}}
    if [[ ${{#reads_right[@]}} -gt 0 ]]; then
        mkdir -p $(dirname {snakemake.output.out}/${{output_right[$i]}})
        mv $TMPDIR/trimmed/${{i}}_R2.fastq.gz {snakemake.output.out}/${{output_right[$i]}}
    fi
done

# AGeNT's trimmer writes no report and keeps no rejected reads
mkdir -p {snakemake.output.out} {snakemake.output.report} {snakemake.output.rejected}
"""
)

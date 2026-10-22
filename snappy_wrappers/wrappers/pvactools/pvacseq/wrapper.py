import csv
import os

from snappy_wrappers.pvactools import bindings_for_container
from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

(input_bindings, input_fns) = bindings_for_container(
    snakemake.input, base_container_path="/inputs", exclude_bind=args["exclude_bind"]
)
(output_bindings, output_fns) = bindings_for_container(
    snakemake.output, base_container_path="/outputs", access="rw", exclude_bind=args["exclude_bind"]
)

if phased := getattr(snakemake.input, "phased", ""):
    phased = f"--phased-proximal-variants-vcf {input_fns['phased']}"

if peptides := getattr(snakemake.input, "peptides", ""):
    peptides = f"--run-reference-proteome-similarity --peptide-fasta {input_fns['peptides']}"
if genes := getattr(snakemake.input, "genes", ""):
    genes = f"--genes-of-interest-file {input_fns['genes']}"

normal = args.get("normal_sample", None)
if normal:
    normal = f"--normal-sample-name {normal}"
else:
    normal = ""

alleles = []
with open(snakemake.input.alleles, "rt") as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        alleles.append(row["HLA Allele"])
alleles = ",".join(sorted(list(set(alleles))))

if "all_class_i" in args["algorithms"] and "all_class_ii" in args["algorithms"]:
    algorithms = "all"
else:
    algorithms = " ".join(args["algorithms"])
class_i_epitope_length = ",".join(map(str, args["lengths"]["class_i"]))
class_ii_epitope_length = ",".join(map(str, args["lengths"]["class_ii"]))

# The epitope sequences that NetChop reads; empty when the tool writes none
fasta = " ".join(str(path) for key, path in snakemake.output.items() if key.startswith("fasta"))

ShellWrapper(snakemake).run(
    r"""
export LC_ALL=C.UTF-8

# snakemake.output.done is needed because pVACtools modules create MHC_class_iI? folders
# so the path to bind is work/mapper.caller.annotator(.filtered)?.tumor_dna/out
# rather than the dirname of snakemake.output.all (which is the MHC_class_II? subdir of the
# bound directory)
out=$(dirname {snakemake.output.done})
rm -rf $out
mkdir -p $out
out=$(realpath $out)

tmp=$(dirname $out)/tmp/pvacseq
rm -rf $tmp
mkdir -p $tmp

scripts=$(dirname $out)/scripts
mkdir -p $scripts

cat << __EOF > $scripts/run_pVACseq.sh
set -x
export TMPDIR=/short_tmp
pvacseq run --n-threads {snakemake.threads} \\
    {normal} \\
    --iedb-install-directory /opt/iedb \\
    --class-i-epitope-length {class_i_epitope_length} --class-ii-epitope-length {class_ii_epitope_length} \\
    {extra_arguments} \\
    {peptides} {genes} {phased} \\
    {input_fns[vcf]} \\
    {args[tumor_sample]} {alleles} {algorithms} \\
    $(dirname {output_fns[done]})
if [[ $? -eq 0 ]]
then
    echo "pvacseq finished without errors, cleaning up tmp"
    rm -rf /short_tmp/*
else
    echo "ERROR during pvacseq"
fi
__EOF
chmod +x $scripts/run_pVACseq.sh

# Short path tmp to avoid AF_UNIX path too long error
apptainer exec \
    --no-home --bind $tmp:/short_tmp:rw --bind $scripts:/scripts:ro \
    {input_bindings} {output_bindings} {snakemake.input[container]} bash /scripts/run_pVACseq.sh

link_missing() {{
    combined=$1
    cls_i=$2
    cls_ii=$3

    if [[ ! -e $combined ]]
    then
        mkdir -p $(dirname $combined)
        if [[ ! -e $cls_i ]]
        then
            mkdir -p $(dirname $cls_i)
            if [[ ! -e $cls_ii ]]
            then
                mkdir -p $(dirname $cls_ii)
                echo "Missing $cls_i and $cls_ii, aborting"
                exit 1
            fi
            touch $cls_i
            ln -sr $cls_ii $combined
        else
            mkdir -p $(dirname $cls_ii)
            touch $cls_ii
            ln -sr $cls_i $combined
        fi
    fi
}}


link_missing {snakemake.output[all.Combined]} {snakemake.output[all.MHC_I]} {snakemake.output[all.MHC_II]}
link_missing {snakemake.output[aggregated.Combined]} {snakemake.output[aggregated.MHC_I]} {snakemake.output[aggregated.MHC_II]}
link_missing {snakemake.output[filtered.Combined]} {snakemake.output[filtered.MHC_I]} {snakemake.output[filtered.MHC_II]}

for path in {fasta}; do
    [[ -e $path ]] || {{ mkdir -p $(dirname $path); touch $path; }}
done

ls -altrR $(dirname {snakemake.output.done})
touch {snakemake.output.done}
"""
)

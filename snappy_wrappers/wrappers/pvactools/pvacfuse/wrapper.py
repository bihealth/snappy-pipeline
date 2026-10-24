import csv

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

args = getattr(snakemake.params, "args", {})

if peptides := getattr(snakemake.input, "peptides", ""):
    peptides = f"--run-reference-proteome-similarity --peptide-fasta {peptides}"
if genes := getattr(snakemake.input, "genes", ""):
    genes = f"--genes-of-interest-file {genes}"

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

# pVACfuse writes the MHC_Class_I, MHC_Class_II and combined directories below $out
out=$(dirname {snakemake.output.done})
rm -rf $out
mkdir -p $out

pvacfuse run --n-threads {snakemake.threads} \
    --iedb-install-directory /opt/iedb \
    --class-i-epitope-length {class_i_epitope_length} --class-ii-epitope-length {class_ii_epitope_length} \
    {extra_arguments} \
    {peptides} {genes} \
    {snakemake.input.fusions} \
    {args[tumor_sample]} {alleles} {algorithms} \
    $out

# Remove lots of tmp file left over by pvacfuse
remove_tmp() {{
    path=$1
    lengths=$2

    lengths=$(echo "$lengths" | tr ',' '\n')

    for l in $lengths; do
        if [[ -d $path/$l/tmp ]]; then
            echo "Clean up $path/$l/tmp"
            rm -rf $path/$l/tmp
        else
            echo "$path/$l/tmp not found"
        fi
    done

    echo "Cleanup complete"
}}
remove_tmp $(dirname {snakemake.output[all.MHC_I]}) "{class_i_epitope_length}"
remove_tmp $(dirname {snakemake.output[all.MHC_II]}) "{class_ii_epitope_length}"

link_missing() {{
    combined=$1
    cls_i=$2
    cls_ii=$3

    if [[ ! -e $combined ]]; then
        mkdir -p $(dirname $combined)
        if [[ ! -e $cls_i ]]; then
            mkdir -p $(dirname $cls_i)
            if [[ ! -e $cls_ii ]]; then
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
    [[ -e $path ]] || {{
        mkdir -p $(dirname $path)
        touch $path
    }}
done

ls -altrR $(dirname {snakemake.output.done})
touch {snakemake.output.done}
"""
)

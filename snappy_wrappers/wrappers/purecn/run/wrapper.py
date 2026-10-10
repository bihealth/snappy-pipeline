# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing CNV using PureCN"""

import os

from snappy_wrappers.snappy_wrapper import ShellWrapper

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

args = getattr(snakemake.params, "args", {})
config = args["config"]

# WARNING- these extra commands cannot contain file paths
extra_commands = " ".join(
    [
        " --{}={}".format(k, v) if v else " --{}".format(k)
        for k, v in config["extra_commands"].items()
    ]
)

# List files that must be accessible from the container
files_to_bind = {
    "vcf": snakemake.input.vcf,
    "mapping_bias": snakemake.input.mapping_bias,
    "normaldb": snakemake.input.normaldb,
    "intervals": snakemake.input.intervals,
}
if "segments" in snakemake.input.keys() and snakemake.input.segments:
    files_to_bind["seg-file"] = snakemake.input.segments
if "log2" in snakemake.input.keys() and snakemake.input.log2:
    files_to_bind["log-ratio-file"] = snakemake.input.log2

# TODO: Put the following in a function (decide where...)
# Replace with full absolute paths
files_to_bind = {k: os.path.realpath(v) for k, v in files_to_bind.items()}
# Directories that mut be bound
dirs_to_bind = {k: os.path.dirname(v) for k, v in files_to_bind.items()}
# List of unique directories to bind: on cluster: <directory> -> from container: /bindings/d<i>)
bound_dirs = {e[1]: e[0] for e in enumerate(list(set(dirs_to_bind.values())))}
# Binding command
bindings = " ".join(["-B {}:/bindings/d{}:ro".format(k, v) for k, v in bound_dirs.items()])
# Path to files from the container
bound_files = {
    k: "/bindings/d{}/{}".format(bound_dirs[dirs_to_bind[k]], os.path.basename(v))
    for k, v in files_to_bind.items()
}

if "seg-file" in bound_files.keys():
    extra_commands += " --seg-file={}".format(bound_files["seg-file"])
if "log-ratio-file" in bound_files.keys():
    extra_commands += " --log-ratio-file={}".format(bound_files["log-ratio-file"])

ShellWrapper(snakemake).run(
    r"""
# Rename PureCN files to snappy conventions
rename() {{
    to=$1
    from=$2
    d=$(dirname $to)
    test -e $d/{args[library_name]}$from
    [ "$d/{args[library_name]}$from" = "$to" ] || mv $d/{args[library_name]}$from $to
}}

outdir=$(dirname {snakemake.output.segments})
mkdir -p $outdir

# Run PureCN with a panel of normals
cmd="/usr/local/bin/Rscript /opt/PureCN/PureCN.R \
    --sampleid {args[library_name]} \
    --tumor {snakemake.input.tumor} \
    --vcf {bound_files[vcf]} \
    --mapping-bias-file {bound_files[mapping_bias]} \
    --normaldb {bound_files[normaldb]} \
    --intervals {bound_files[intervals]} \
    --genome {config[genome_name]} \
    --out $outdir --out-vcf --force \
    --seed {config[seed]} --parallel --cores {snakemake.threads} \
    {extra_commands}
"
apptainer exec --home $PWD {bindings} {snakemake.input.container} $cmd

rename {snakemake.output.segments} _dnacopy.seg
rename {snakemake.output.ploidy} .csv
rename {snakemake.output.pvalues} _amplification_pvalues.csv
rename {snakemake.output.vcf} .vcf.gz
rename {snakemake.output.vcf_tbi} .vcf.gz.tbi
rename {snakemake.output.loh} _loh.csv

# Fix chromosome names (https://github.com/lima1/PureCN/issues/331)
vcf_chrnames=$(zgrep '^##contig=<ID=' {snakemake.input.vcf} | sed -re "s/^##contig=<ID=([^,]*),.*/\1/" | sort | uniq | grep -E "^(chr)?([0-9]+|[XY])$")
n_prefix=0
n_tot=0
for chrname in $vcf_chrnames
do
    ((n_tot=n_tot + 1))
    n=$(echo "$chrname" | grep -c "^$chrname$" || true)
    n_prefix=$((n_prefix + $n))
done
[ $n_prefix -eq 0 ] || [ $n_prefix -eq $n_tot ]

if [[ $n_prefix -eq $n_tot ]]
then
    pgm='{{
        chr=$2
        if(chr=="23"){{chr="X"}}
        if(chr=="24"){{chr="Y"}}
        if($0 ~ /^ID\t/){{
            print $0
        }} else {{
            printf"%s\tchr%s\t%d\t%d\t%d\t%f\t%d\n",$1,chr,int($3+0.5),int($4+0.5),int($5+0.5),$6,int($7+0.5)
        }}
    }}'
else
    pgm='{{
        chr=$2
        if(chr=="23"){{chr="X"}}
        if(chr=="24"){{chr="Y"}}
        if($0 ~ /^ID\t/){{
            print $0
        }} else {{
            printf"%s\t%s\t%d\t%d\t%d\t%f\t%d\n",$1,chr,int($3+0.5),int($4+0.5),int($5+0.5),$6,int($7+0.5)
        }}
    }}'
fi
mv {snakemake.output.segments} $TMPDIR/segments.seg
awk -F'\t' "$pgm" $TMPDIR/segments.seg > {snakemake.output.segments}
"""
)

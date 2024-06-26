# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for computing CNV using PureCN"""

import os

from snakemake import shell

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

step = snakemake.config["pipeline_step"]["name"]
config = snakemake.config["step_config"][step]["purecn"]

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
    "mapping_bias": config["path_mapping_bias"],
    "normaldb": config["path_panel_of_normals"],
    "intervals": config["path_intervals"],
}
if "snp_blacklist" in config.keys() and config["snp_blacklist"]:
    files_to_bind["snp-blacklist"] = config["snp_blacklist"]
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

if "snp-blacklist" in bound_files.keys():
    extra_commands += " --snp-blacklist={}".format(bound_files["snp-blacklist"])
if "seg-file" in bound_files.keys():
    extra_commands += " --seg-file={}".format(bound_files["seg-file"])
if "log-ratio-file" in bound_files.keys():
    extra_commands += " --log-ratio-file={}".format(bound_files["log-ratio-file"])

shell.executable("/bin/bash")

shell(
    r"""
set -x

# Also pipe everything to log file
if [[ -n "{snakemake.log.log}" ]]; then
    if [[ "$(set +e; tty; set -e)" != "" ]]; then
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        exec &> >(tee -a "{snakemake.log.log}" >&2)
    else
        rm -f "{snakemake.log.log}" && mkdir -p $(dirname {snakemake.log.log})
        echo "No tty, logging disabled" >"{snakemake.log.log}"
    fi
fi

# Write out information about conda installation.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}

# Setup auto-cleaned tmpdir
export tmpdir=$(mktemp -d)
trap "rm -rf $tmpdir" EXIT

# Compute md5 checksum
md5() {{
    fn=$1
    d=$(dirname $fn)
    f=$(basename $fn)
    pushd $d 1> /dev/null 2>&1
    checksum=$(md5sum $f)
    popd 1> /dev/null 2>&1
    echo "$checksum"
}}

# Rename PureCN files to snappy conventions
rename() {{
    to=$1
    d=$(dirname $to)
    f=$(basename $to)
    echo $f | grep -q "^{snakemake.wildcards[mapper]}.purecn.{snakemake.wildcards[library_name]}"
    from=$(echo $f | sed -e "s/^{snakemake.wildcards[mapper]}.purecn.//")
    test -e $d/$from
    mv $d/$from $to
}}

outdir=$(dirname {snakemake.output.segments})
mkdir -p $outdir

# Run PureCN with a panel of normals
cmd="/usr/local/bin/Rscript /opt/PureCN/PureCN.R \
    --sampleid {snakemake.wildcards[library_name]} \
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
apptainer exec --home $PWD {bindings} {config[path_container]} $cmd

rename {snakemake.output.segments}
rename {snakemake.output.ploidy}
rename {snakemake.output.pvalues}
rename {snakemake.output.vcf}
rename {snakemake.output.vcf_tbi}
rename {snakemake.output.loh}

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
mv {snakemake.output.segments} $tmpdir/segments.seg
awk -F'\t' "$pgm" $tmpdir/segments.seg > {snakemake.output.segments}

md5 {snakemake.output.segments} > {snakemake.output.segments_md5}
md5 {snakemake.output.ploidy} > {snakemake.output.ploidy_md5}
md5 {snakemake.output.pvalues} > {snakemake.output.pvalues_md5}
md5 {snakemake.output.vcf} > {snakemake.output.vcf_md5}
md5 {snakemake.output.vcf_tbi} > {snakemake.output.vcf_tbi_md5}
md5 {snakemake.output.loh} > {snakemake.output.loh_md5}
"""
)

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

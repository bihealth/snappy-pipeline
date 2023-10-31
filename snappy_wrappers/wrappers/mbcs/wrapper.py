# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MTB-aware exome data
"""

import os
import shutil
import sys
import tempfile

from snakemake import shell

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import run_snakemake

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

shell.executable("/bin/bash")

# Rules templates =============================================================

# Main rule: rule all ---------------------------------------------------------
head_rule = r"""
shell.prefix("set -xeuo pipefail; ")

rule all:
    input:
        "{output_bam}"
"""

# Trim UMIs & map one pair of fastq files -------------------------------------
paired_rule = r"""
rule trim_{iPair}:
    input:
        r1 = "{r1}",
        r2 = "{r2}"
    output:
        r1 = "trimmed/{name}_R1.fastq.gz",
        r2 = "trimmed/{name}_R2.fastq.gz"
    log:
        "log/{name}_trim.log"
    resources:
        mem_mb = "64000M",
        time = "24:00:00"
    shell:
        r'''
            java -Xmx60G -jar {trimmer} \
                -{lib_prep_type} {extra_args} \
                -out "trimmed/{name}" \
                -fq1 {{input.r1}} -fq2 {{input.r2}}
        '''

rule map_{iPair}:
    input:
        r1 = "trimmed/{name}_R1.fastq.gz",
        r2 = "trimmed/{name}_R2.fastq.gz"
    output:
        bam = "mapped/{name}.bam"
    log:
        "log/{name}_map.log"
    resources:
        mem_mb = "64000M",
        time = "24:00:00"
    threads: 16
    params:
        indices = "{indices}"
    shell:
        r'''
            seqtk mergepe <(zcat {{input.r1}}) <(zcat {{input.r2}}) \
                | {mapper} mem {{params.indices}} -t {{threads}} -p -R '@RG\tID:{name}\tSM:{sample}\tPL:ILLUMINA' -C /dev/stdin \
                | samtools sort -O BAM -o {{output.bam}}
            samtools index {{output.bam}}
        '''
"""

# Trim & map one fastq file (unpaired) ----------------------------------------
unpaired_rule = r"""
rule trim_{iPair}:
    input:
        r1 = "{r1}"
    output:
        r1 = "trimmed/{name}_R1.fastq.gz"
    log:
        "log/{name}_trim.log"
    resources:
        mem_mb = "64000M",
        time = "24:00:00"
    shell:
        r'''
            java -Xmx60G -jar {trimmer} -fq1 {{input.r1}}
        '''

rule map_{iPair}:
    input:
        r1 = "trimmed/{name}_R1.fastq.gz"
    output:
        bam = "mapped/{name}.bam"
    log:
        "log/{name}_map.log"
    resources:
        mem_mb = "64000M",
        time = "24:00:00"
    threads: {threads}
    params:
        indices = "{indices}"
    shell:
        r'''
            seqtk mergepe <(zcat {{input.r1}}) \
                | {mapper} mem {{params.indices}} -t {{threads}} -p -R '@RG\tID:{name}\tSM:{sample}\tPL:ILLUMINA' -C /dev/stdin \
                | samtools sort -O BAM -o {{output.bam}}
            samtools index {{output.bam}}
        '''
"""

# Merge, mark duplicated and base quality recalibration on all data -----------
merge_rule = r"""
rule merge:
    input:
        pairs = [{pairs}]
    output:
        bam = "merged.bam"
    log:
        "log/merge.log"
    shell:
        r'''
            samtools merge -O BAM -o {{output.bam}} -f -s 1234567 {{input.pairs}}
            samtools index {{output.bam}}
        '''

rule mark:
    input:
        bam = "merged.bam"
    output:
        bam = "marked.bam"
    log:
        "log/marked.log"
    resources:
        mem_mb = "64000M",
        time = "24:00:00"
    shell:
        r'''
            java -Xmx60G -jar {marker} \
                {extra_args} -c={consensus_mode} \
                -o={{output.bam}} \
                -f {input_filer_args} -F {consensus_filter_args} \
                {{input.bam}}

            samtools index {{output.bam}}
        '''
    
rule bqsr:
    input:
        bam = "marked.bam"
    output:
        tbl = "bqsr.tbl"
    log:
        "log/bqsr.log"
    params:
        reference = "{reference}",
        common_sites = "{common_sites}"
    shell:
        r'''
            gatk BaseRecalibrator \
                -I {{input.bam}} \
                -R {{params.reference}} --known-sites {{params.common_sites}} \
                -O {{output.tbl}}
        '''

rule apply:
    input:
        bam = "marked.bam",
        tbl = "bqsr.tbl"
    output:
        bam = "{output_bam}"
    log:
        "log/apply.log"
    params:
        reference = "{reference}"
    shell:
        r'''
            gatk ApplyBQSR \
                -I {{input.bam}} --bqsr-recal-file {{input.tbl}} \
                -R {{params.reference}} \
                -O {{output.bam}}

            samtools index {{output.bam}}
        '''
"""

# Main script =================================================================

# Input fastqs are passed through snakemake.params.
# snakemake.input is a .done file touched after linking files in.
input_left = snakemake.params.args["input"]["reads_left"]
input_right = snakemake.params.args["input"].get("reads_right", "")

config = snakemake.config["step_config"]["ngs_mapping"]["mbcs"]

# Group mate pairs
r1s = input_left.copy()

pairs = {}
for r2 in input_right:
    f = os.path.basename(r2).replace("_R2", "_R1")
    r1 = None
    for x in r1s:
        if os.path.basename(x) == f:
            r1 = x
            break
    assert r1 is not None
    name = f[: f.index("_R1")]
    assert name not in pairs.keys()
    pairs[name] = (os.path.abspath(r1), os.path.abspath(r2))
    r1s.remove(r1)

for r1 in r1s:
    f = os.path.basename(r1)
    assert "_R1" in f
    name = f[: f.index("_R1")]
    assert name not in pairs.keys()
    pairs[name] = (os.path.abspath(r1),)

# Create Snakefile
snakefile = []

rule = head_rule.format(output_bam=os.path.abspath(snakemake.output.bam))
snakefile.append(rule)

mapper = config["mapping_tool"]
tool = config["mbc_tool"]
for iPair, name in enumerate(pairs.keys()):
    kwargs = {
        "iPair": iPair,
        "name": name,
        "r1": pairs[name][0],
        "r2": pairs[name][1],
        "trimmer": config[tool]["prepare"]["path"],
        "lib_prep_type": config[tool]["prepare"]["lib_prep_type"],
        "extra_args": " ".join(config[tool]["prepare"]["extra_args"]),
        "mapper": mapper if mapper != "bwa_mem2" else "bwa-mem2",
        "indices": snakemake.config["step_config"]["ngs_mapping"][mapper]["path_index"],
        "sample": snakemake.params.args["sample_name"],
        "threads": snakemake.config["step_config"]["ngs_mapping"][mapper]["num_threads_align"],
    }
    if len(pairs[name]) == 2:
        rule = paired_rule.format(**kwargs)
    else:
        rule = unpaired_rule.format(**kwargs)
    snakefile.append(rule)

rule = merge_rule.format(
    pairs='"' + '", "'.join(["mapped/{}.bam".format(name) for name in pairs.keys()]) + '"',
    marker=config[tool]["mark_duplicates"]["path"],
    consensus_mode=config[tool]["mark_duplicates"]["consensus_mode"],
    input_filer_args=" ".join(config[tool]["mark_duplicates"]["input_filter_args"]),
    consensus_filter_args=" ".join(config[tool]["mark_duplicates"]["consensus_filter_args"]),
    extra_args=" ".join(config[tool]["mark_duplicates"]["extra_args"]),
    reference=snakemake.config["static_data_config"]["reference"]["path"],
    common_sites=snakemake.config["step_config"]["ngs_mapping"]["mbcs"]["agent"]["common_variants"],
    output_bam=os.path.abspath(snakemake.output.bam),
)
snakefile.append(rule)

snakefile = "\n\n# =======================================================\n\n".join(snakefile)

# Create temp directory, go there, write Snakefile and run
pwd = os.getcwd()

# import pdb; pdb.set_trace()

tempdir = tempfile.mkdtemp()
os.chdir(tempdir)

# os.mkdir(os.path.join(tempdir, "slurm_log"), mode=0o750)
with open("Snakefile", "wt") as f:
    f.write(snakefile)

slurm_config = {
    "num_jobs": 10,
    "max_jobs_per_second": 10,
    "max_status_checks_per_second": 1,
    "job_name_token": "ngs_mapping_mbcs",
    "profile": os.getenv("SNAPPY_PIPELINE_SNAKEMAKE_PROFILE"),
}
config["use_profile"] = True
config["restart_times"] = 1
run_snakemake(
    config,
    **slurm_config,
)

os.chdir(pwd)

# Finish: write stats, logs & links
shell(
    r"""
set -x

# Write out information about conda and save a copy of the wrapper with picked variables
# as well as the environment.yaml file.
conda list >{snakemake.log.conda_list}
conda info >{snakemake.log.conda_info}
md5sum {snakemake.log.conda_list} >{snakemake.log.conda_list_md5}
md5sum {snakemake.log.conda_info} >{snakemake.log.conda_info_md5}
cp {__real_file__} {snakemake.log.wrapper}
md5sum {snakemake.log.wrapper} >{snakemake.log.wrapper_md5}
cp $(dirname {__file__})/environment.yaml {snakemake.log.env_yaml}
md5sum {snakemake.log.env_yaml} >{snakemake.log.env_yaml_md5}

# Concatenate all log files into main log
for fn in $(ls -tr {tempdir}/slurm_log/*/*.log)
do
    f=$(basename $fn)
    echo "# ==============================================================" >> {snakemake.log.log}
    echo "# Log from \"$f\"" >> {snakemake.log.log}
    echo "# ==============================================================" >> {snakemake.log.log}
    cat $fn >> {snakemake.log.log}
done
md5sum {snakemake.log.log} > {snakemake.log.log_md5}

# Copy BQSR table just in case (not on output)
mv {tempdir}/bqsr.tbl $(dirname {snakemake.output.bam})

# Build MD5 files
pushd $(dirname {snakemake.output.bam})
md5sum $(basename {snakemake.output.bam}) > $(basename {snakemake.output.bam}).md5
md5sum $(basename {snakemake.output.bam_bai}) > $(basename {snakemake.output.bam_bai}).md5
popd

# QC Report ---------------------------------------------------------------------------------------

# gather statistics from BAM file
# TODO: use pipes for only reading once from disk?
samtools stats    {snakemake.output.bam} > {snakemake.output.report_bamstats_txt}
samtools flagstat {snakemake.output.bam} > {snakemake.output.report_flagstats_txt}
samtools idxstats {snakemake.output.bam} > {snakemake.output.report_idxstats_txt}

# Build MD5 files for the reports
md5sum {snakemake.output.report_bamstats_txt} > {snakemake.output.report_bamstats_txt_md5}
md5sum {snakemake.output.report_flagstats_txt} >{snakemake.output.report_flagstats_txt_md5}
md5sum {snakemake.output.report_idxstats_txt} > {snakemake.output.report_idxstats_txt_md5}

# Create output links -----------------------------------------------------------------------------

for path in {snakemake.output.output_links}; do
  dst=$path
  src=work/${{dst#output/}}
  ln -sr $src $dst
done
"""
)

shutil.rmtree(tempdir)

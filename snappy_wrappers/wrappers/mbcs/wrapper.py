# -*- coding: utf-8 -*-
"""CUBI+Snakemake wrapper code for MTB-aware exome data"""

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

args = getattr(snakemake.params, "args", {})

# Helper functions ------------------------------------------------------------
def pair_fastq_files(input_left, input_right):
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

    return pairs


# Read snakemake input --------------------------------------------------------
input_left = args["input"]["reads_left"]
input_right = args["input"].get("reads_right", "")

config = args["config"]
mapper = config["mapping_tool"]
mapper_config = args["mapper_config"]
if mapper == "bwa_mem2":
    mapper = "bwa-mem2"
if config["use_barcodes"]:
    barcoder = config["barcode_tool"]
    config_barcodes = args["barcode_config"]
if config["recalibrate"]:
    config_bqsr = args["bqsr_config"]

# Group fastq files by lane ---------------------------------------------------
pairs = pair_fastq_files(input_left, input_right)

# Rules templates -------------------------------------------------------------
head_rule = r"""
shell.prefix("set -xeuo pipefail; ")

rule all:
    input:
        "{output_bam}"
"""

generic_rule = r"""
rule {rule}:
    input:
        {input}
    output:
        {output}
    threads: {threads}
    resources:
        mem_mb = "{mem}",
        time = "{time}",
        partition = "medium"
    shell:
        r'''
{cmds}
        '''
"""

# Format input for each lane (r1 = ... & optionally r2 = ...) -----------------
inputs = {}
for iPair, name in enumerate(pairs.keys()):
    pair = pairs[name]
    in_ = 'r1 = "{}"'.format(pair[0])
    if len(pair) == 2:
        in_ += ',\n        r2 = "{}"'.format(pair[1])
    inputs[name] = in_

# Create snakefile chunks =====================================================
snakefile = [head_rule.format(output_bam=os.path.abspath(snakemake.output.bam))]

# Trimming barcodes -----------------------------------------------------------
if config["use_barcodes"]:
    for iPair, name in enumerate(pairs.keys()):
        pair = pairs[name]

        out = 'r1 = "trimmed/{}_R1.fastq.gz"'.format(name)
        fq2 = ""
        if len(pair) == 2:
            out += ',\n        r2 = "trimmed/{}_R2.fastq.gz"'.format(name)
            fq2 = "-fq2 {input.r2}"

        cmd = (
            "            java -Xmx60G -jar {trimmer} \\\n"
            "                -{lib_prep_type} {extra_args} \\\n"
            '                -out "trimmed/{name}" \\\n'
            "                -fq1 {{input.r1}} {fq2}"
        ).format(
            trimmer=config_barcodes["prepare"]["path"],
            lib_prep_type=config_barcodes["prepare"]["lib_prep_type"],
            extra_args=" ".join(config_barcodes["prepare"]["extra_args"]),
            name=name,
            fq2=fq2,
        )

        kwargs = {
            "rule": f"trim_{iPair}",
            "input": inputs[name],
            "output": out,
            "cmds": cmd,
            "mem": "64000M",
            "time": "24:00:00",
            "threads": "1",
        }
        snakefile.append(generic_rule.format(**kwargs))

        inputs[name] = out

# Mapping chunk ===============================================================

# Prepare mapping command (many steps, with pipes) ----------------------------

# Input depends on presence of R2
cmd = "            {zcat}"

# Optionally trim adapters
if mapper_config["trim_adapters"]:
    cmd += "            | trimadap-mt -p {} \\\n".format(mapper_config["num_threads_trimming"])

# Map with read group (bwa & bwa-mem2 only)
cmd += (
    "            | {mapper} mem \\\n"
    "                {indices} \\\n"
    "                -R '@RG\tID:{{name}}\tSM:{sample_name}\tPL:ILLUMINA' \\\n"
    "                -p {extra_args} -t {threads} /dev/stdin \\\n"
).format(
    mapper=mapper,
    indices=mapper_config["path_index"],
    sample_name=args["sample_name"],
    extra_args=" ".join(mapper_config.get("extra_args", [])),
    threads=mapper_config["num_threads_align"],
)

# Optionally mask duplicates
if mapper_config["mask_duplicates"]:
    cmd += "            | samblaster --addMateTags \\\n"

# Sort output bam - end of the pipe
cmd += (
    "            | samtools sort -T {{tmp}} \\\n"
    "                -m {memory} -@ {threads} -O BAM -o {{{{output.bam}}}} /dev/stdin\n"
).format(
    memory=mapper_config["memory_bam_sort"],
    threads=mapper_config["num_threads_bam_sort"],
)

# Index output bam
cmd += "            samtools index {{output.bam}}"

# Loop over lanes -------------------------------------------------------------
for iPair, name in enumerate(pairs.keys()):
    pair = pairs[name]

    if len(pair) == 2:
        zcat = "seqtk mergepe <(zcat {input.r1}) <(zcat {input.r2}) \\\n"
    else:
        zcat = "zcat {input.r1} \\\n"

    if len(pairs.keys()) > 1 or config["use_barcodes"] or config["recalibrate"]:
        out = f"mapped/{name}.bam"
    else:
        out = os.path.abspath(snakemake.output.bam)

    kwargs = {
        "rule": f"map_{iPair}",
        "input": inputs[name],
        "output": f'bam = "{out}"',
        "cmds": cmd.format(zcat=zcat, name=name, out=out, tmp=f"sort/{name}"),
        "threads": mapper_config["num_threads_align"] + mapper_config["num_threads_bam_sort"],
        "mem": "64000M",
        "time": "24:00:00",
    }
    snakefile.append(generic_rule.format(**kwargs))

# Merge bam files when more than one lane -------------------------------------
if len(pairs.keys()) > 1:
    if config["use_barcodes"] or config["recalibrate"]:
        out = "merged.bam"
    else:
        out = os.path.abspath(snakemake.output.bam)
    cmd = (
        "            samtools merge -O BAM -o {output.bam} -f -s 1234567 {input.bams} \n"
        "            samtools index {output.bam}"
    )
    rule = generic_rule.format(
        rule="merge",
        input="bams = [{}]".format(
            ", ".join(['"mapped/{}.bam"'.format(name) for name in pairs.keys()])
        ),
        output=f'bam = "{out}"',
        cmds=cmd,
        threads="1",
        mem="16000M",
        time="4:00:00",
    )
    snakefile.append(rule)
in_ = out

# Mark duplicates with barcodes -----------------------------------------------
if config["use_barcodes"]:
    if config["recalibrate"]:
        out = "marked.bam"
    else:
        out = os.path.abspath(snakemake.output.bam)

    cmd = (
        "            java -Xmx60G -jar {marker} {extra_args} \\\n"
        "                --bed-file={baits} \\\n"
        "                -c={consensus_mode} \\\n"
        "                -o={{output.bam}} \\\n"
        "                -f {input_filter_args} -F {consensus_filter_args} \\\n"
        "                {{input.bam}} \n"
        "            samtools index {{output.bam}}"
    ).format(
        marker=config_barcodes["mark_duplicates"]["path"],
        consensus_mode=config_barcodes["mark_duplicates"]["consensus_mode"],
        baits=config_barcodes["mark_duplicates"]["path_baits"],
        extra_args=" ".join(config_barcodes["mark_duplicates"]["extra_args"]),
        input_filter_args=" ".join(config_barcodes["mark_duplicates"]["input_filter_args"]),
        consensus_filter_args=" ".join(config_barcodes["mark_duplicates"]["consensus_filter_args"]),
    )

    kwargs = {
        "rule": "mark",
        "input": f'bam = "{in_}"',
        "output": f'bam = "{out}"',
        "cmds": cmd,
        "mem": "96000M",
        "time": "24:00:00",
        "threads": "1",
    }
    snakefile.append(generic_rule.format(**kwargs))
in_ = out

# Base quality recalibration --------------------------------------------------
if config["recalibrate"]:
    cmd = (
        "            gatk BaseRecalibrator \\\n"
        "                -I {{input.bam}} \\\n"
        "                -R {reference} --known-sites {common_sites} \\\n"
        "                -O {{output.tbl}}"
    ).format(
        reference=args["reference"],
        common_sites=config_bqsr["common_variants"],
    )

    kwargs = {
        "rule": "bqsr",
        "input": f'bam = "{in_}"',
        "output": 'tbl = "bqsr.tbl"',
        "cmds": cmd,
        "mem": "32000M",
        "time": "24:00:00",
        "threads": "1",
    }
    snakefile.append(generic_rule.format(**kwargs))

    cmd = (
        "            gatk ApplyBQSR \\\n"
        "                -I {{input.bam}} --bqsr-recal-file {{input.tbl}} \\\n"
        "                -R {reference} \\\n"
        "                -O {{output.bam}} \n"
        "            samtools index {{output.bam}}"
    ).format(reference=args["reference"])

    kwargs = {
        "rule": "apply",
        "input": f'bam = "{in_}",\n        tbl = "bqsr.tbl"',
        "output": 'bam = "{}"'.format(os.path.abspath(snakemake.output.bam)),
        "cmds": cmd,
        "mem": "32000M",
        "time": "24:00:00",
        "threads": "1",
    }
    snakefile.append(generic_rule.format(**kwargs))

# Create snakefile from chunks
snakefile = "\n\n# =======================================================\n\n".join(snakefile)

# Start snakemake in temporary directory ======================================
pwd = os.getcwd()

tempdir = tempfile.mkdtemp()
os.chdir(tempdir)
os.mkdir("sort", mode=0o750)

with open("Snakefile", "wt") as f:
    f.write(snakefile)

slurm_config = {
    "cores": mapper_config["num_threads_align"] + mapper_config["num_threads_bam_sort"],
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

# Finish: write stats, logs & links -------------------------------------------
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
jobid=$(ls {tempdir}/slurm_log)
fns=$(ls -tr {tempdir}/slurm_log/$jobid/*.log)
ls -al $fns
for fn in $fns
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

rm -rf {tempdir}
"""
)

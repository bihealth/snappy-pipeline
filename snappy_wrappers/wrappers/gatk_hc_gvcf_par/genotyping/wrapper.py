import json
import os
import sys
import tempfile
import textwrap

base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.tools.genome_windows import yield_regions
from snappy_wrappers.wrapper_parallel import (
    ResourceUsage,
    SgeResourceUsageConverter,
    gib,
    hours,
    in_working_dir,
    run_snakemake,
)

# The step and caller key to use.
step_key = snakemake.params.step_key
caller_key = snakemake.params.caller_key

# TODO: call on overlapping windows, on merge make unique

# Naming clash limbo...
snake_job = snakemake
del snakemake
from snakemake import shell
from snakemake.io import Namedlist

snakemake = snake_job

# Fix-Up Paths ------------------------------------------------------------------------------------

snakemake.input = Namedlist(map(os.path.realpath, snakemake.input))
snakemake.output.vcf = os.path.realpath(snakemake.output.vcf)
snakemake.output.vcf_md5 = os.path.realpath(snakemake.output.vcf_md5)
snakemake.output.tbi = os.path.realpath(snakemake.output.tbi)
snakemake.output.tbi_md5 = os.path.realpath(snakemake.output.tbi_md5)

# Perform Splitting -------------------------------------------------------------------------------

# Split chromosome (lengths from FAI file) into windows of configured length
fai_path = snakemake.config["static_data_config"]["reference"]["path"] + ".fai"
window_size = snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"]["window_length"]
ignore_chroms = snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"]["ignore_chroms"]
with open(fai_path, "rt") as fai_file:
    regions = list(yield_regions(fai_file, window_size, ignore_chroms=ignore_chroms))


# Generate Snakefile chunk-by-chunk ---------------------------------------------------------------
result_files = ["job_out.{jobno}.d/.done".format(jobno=jobno) for jobno in range(len(regions))]
chunks = [
    textwrap.dedent(
        r"""
    shell.executable("/bin/bash")
    shell.prefix("set -ex;")

    configfile: 'config.json'

    localrules: all

    rule all:
        input: [{all_results}]
    """
    )
    .lstrip()
    .format(all_results=", ".join(map(repr, result_files)))
]


def esc(s):
    return s.replace("{", "{{").replace("}", "}}")


#: Extensions to generate
key_ext = {
    "vcf": "vcf.gz",
    "vcf_md5": "vcf.gz.md5",
    "tbi": "vcf.gz.tbi",
    "tbi_md5": "vcf.gz.tbi.md5",
}
if any(".combine_gvcf." in x for x in snakemake.input):
    resources = ResourceUsage(cores=2, memory=gib(30.0), duration=hours(7.5))
else:
    resources = ResourceUsage(cores=2, memory=gib(12.0), duration=hours(7.5))
for jobno, region in enumerate(regions):
    params = dict(snakemake.params)
    if "whole_cohort" in str(snakemake.output):
        chrom_name = region.chrom.replace(".", "_")
        input_gvcf = [
            fname
            for fname in snakemake.input
            if "combine_gvcf.{}.".format(chrom_name) in fname and fname.endswith(".vcf.gz")
        ]
    else:
        input_gvcf = snakemake.input
    params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
    output = {
        key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
        for key, ext in key_ext.items()
    }
    chunks.append(
        textwrap.dedent(
            r"""
    rule chunk_{jobno}:
        input:
            {input_gvcf},
        output:
            touch("job_out.{jobno}.d/.done"),
            **{output}
        params:
            **{params}
        wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/gatk_hc_gvcf/genotyping'


    cluster_config['chunk_{jobno}'] = {resources}
    """
        ).format(
            jobno=jobno,
            params=repr(params),
            output=repr(output),
            input_gvcf=repr(input_gvcf),
            wrapper_prefix="file://" + base_dir,
            resources=repr(SgeResourceUsageConverter(resources).to_res_dict()),
        )
    )

tmpdir = tempfile.mkdtemp("snake_par")
with in_working_dir(tmpdir, print_chdir=True):
    with open("Snakefile", "wt") as snakefile:
        print("\n\n".join(chunks), file=sys.stderr)
        print("\n\n".join(chunks), file=snakefile)

# Write out config file
with in_working_dir(tmpdir, print_chdir=True):
    with open("config.json", "wt") as configfile:
        json.dump(snakemake.config, configfile)
    with open("config.json", "rt") as configfile:
        print(configfile.read(), file=sys.stderr)


# Launch execution of Snakefile
with in_working_dir(tmpdir, print_chdir=True):
    if "whole_cohort" in str(snakemake.output):
        num_jobs = snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"].get(
            "num_jobs_genotype_cohort"
        )
    else:
        num_jobs = None  # fall back to pedigree job count
    run_snakemake(
        snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"],
        num_jobs=num_jobs,
        max_jobs_per_second=snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"][
            "max_jobs_per_second"
        ],
        max_status_checks_per_second=snakemake.config["step_config"]["variant_calling"][
            "gatk_hc_gvcf"
        ]["max_status_checks_per_second"],
        job_name_token="gatk_hc_gvcf_par_genotyping",
        drmaa_snippet=(
            snakemake.config["step_config"][step_key][caller_key]["drmaa_snippet"]
            or snakemake.config["step_config"][step_key]["drmaa_snippet"]
        ),
    )

# Join results
max_jobno = len(regions) - 1
with in_working_dir(tmpdir, print_chdir=True):
    shell(
        textwrap.dedent(
            r"""
    set -euo pipefail

    # Also pipe everything to log file
    if [[ -n "{snakemake.log}" ]]; then
        if [[ "$(set +e; tty; set -e)" != "" ]]; then
            rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
            exec &> >(tee -a "{snakemake.log}" >&2)
        else
            rm -f "{snakemake.log}" && mkdir -p $(dirname {snakemake.log})
            echo "No tty, logging disabled" >"{snakemake.log}"
        fi
    fi

    outdir=$(basename {snakemake.output.vcf})

    mkdir -p output

    # take first header -------------------------------------------------------
    set +e
    zgrep '^#' job_out.0.d/out/tmp_0.vcf.gz > output/out.vcf
    set -e

    # concatenate files -------------------------------------------------------
    bcftools concat \
        -o output/out.vcf.gz \
        -O z \
        $(for jobno in {{0..{max_jobno}}}; do echo job_out.$jobno.d/out/tmp_$jobno.vcf.gz; done)

    tabix -f output/out.vcf.gz

    pushd output && for f in *; do md5sum $f >$f.md5; done && popd

    # move to output directory ------------------------------------------------
    mkdir -p $(dirname {snakemake.output.vcf})
    mv output/out.vcf.gz {snakemake.output.vcf}
    mv output/out.vcf.gz.md5 {snakemake.output.vcf_md5}
    mv output/out.vcf.gz.tbi {snakemake.output.tbi}
    mv output/out.vcf.gz.tbi.md5 {snakemake.output.tbi_md5}
    """
        ).lstrip()
    )

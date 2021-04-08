# isort:skip_file
import json
import os
import sys
import tempfile
import textwrap

base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    ResourceUsage,
    SgeResourceUsageConverter,
    gib,
    hours,
    in_working_dir,
    run_snakemake,
)

# TODO: call on overlapping windows, on merge make unique

# Naming clash limbo...
snake_job = snakemake
del snakemake
from snakemake.io import Namedlist

snakemake = snake_job

# TODO: make calling workflow expect per-contig gvcf files

# Fix-Up Paths ------------------------------------------------------------------------------------

snakemake.input = Namedlist(map(os.path.realpath, snakemake.input))
snakemake.output.vcf = Namedlist(map(os.path.realpath, snakemake.output.vcf))
snakemake.output.vcf_md5 = Namedlist(map(os.path.realpath, snakemake.output.vcf_md5))
snakemake.output.tbi = Namedlist(map(os.path.realpath, snakemake.output.tbi))
snakemake.output.tbi_md5 = Namedlist(map(os.path.realpath, snakemake.output.tbi_md5))

# Generate Snakefile chunk-by-chunk ---------------------------------------------------------------

regions = snakemake.params.args["genome_regions"]

result_files = [
    "job_out.{chrom}.{region_no}.d/.done".format(chrom=chrom.replace(".", "_"), region_no=region_no)
    for chrom, chrom_regions in regions.items()
    for region_no, _ in enumerate(chrom_regions)
]
chunks_work = [
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
    "vcf": "g.vcf.gz",
    "vcf_md5": "g.vcf.gz.md5",
    "tbi": "g.vcf.gz.tbi",
    "tbi_md5": "g.vcf.gz.tbi.md5",
}
resources = ResourceUsage(cores=2, memory=gib(14.0), duration=hours(7.5))
for chrom, chrom_regions in regions.items():
    chrom_name = chrom.replace(".", "_")
    for region_no, region in enumerate(chrom_regions):
        params = dict(snakemake.params)
        region_str = "{}:{:,}-{:,}".format(region["chrom"], region["begin"] + 1, region["end"])
        params.setdefault("args", {}).update({"intervals": [region_str]})
        output = {
            key: "job_out.{chrom}.{region_no}.d/out/tmp_{chrom}.{region_no}.{ext}".format(
                chrom=chrom_name, region_no=region_no, ext=ext
            )
            for key, ext in key_ext.items()
        }
        chunks_work.append(
            textwrap.dedent(
                r"""
        rule chunk_{chrom}_{region_no}:
            input:
                {input_gvcf},
            output:
                touch('job_out.{chrom}.{region_no}.d/.done'),
                **{output}
            params:
                **{params}
            wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/gatk_hc_gvcf/combine_gvcf'


        cluster_config['chunk_{chrom}_{region_no}'] = {resources}
        """
            ).format(
                chrom=chrom_name,
                region_no=region_no,
                params=repr(params),
                output=repr(output),
                input_gvcf=repr(list(snakemake.input)),
                wrapper_prefix="file://" + base_dir,
                resources=repr(SgeResourceUsageConverter(resources).to_res_dict()),
            )
        )

tmpdir_work = tempfile.mkdtemp("snake_par")
with in_working_dir(tmpdir_work, print_chdir=True):
    with open("Snakefile", "wt") as snakefile:
        print("\n\n".join(chunks_work), file=sys.stderr)
        print("\n\n".join(chunks_work), file=snakefile)

# Write out config file
with in_working_dir(tmpdir_work, print_chdir=True):
    with open("config.json", "wt") as configfile:
        json.dump(snakemake.config, configfile)
    with open("config.json", "rt") as configfile:
        print(configfile.read(), file=sys.stderr)


# Launch execution of Snakefile
with in_working_dir(tmpdir_work, print_chdir=True):
    run_snakemake(
        snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"],
        num_jobs=snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"].get(
            "num_jobs_combine_gvcf_cohort"
        ),
        max_jobs_per_second=snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"][
            "max_jobs_per_second"
        ],
        max_status_checks_per_second=snakemake.config["step_config"]["variant_calling"][
            "gatk_hc_gvcf"
        ]["max_status_checks_per_second"],
        job_name_token="gatk_hc_gvcf_par_combine_gvcf_work",
        drmaa_snippet=(
            snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"]["drmaa_snippet"]
            or snakemake.config["step_config"]["variant_calling"]["drmaa_snippet"]
        ),
    )

# Generate Snakefile chunk-by-chunk ---------------------------------------------------------------
# for joining results

result_files_join = [
    "job_out.{chrom}.d/.done".format(chrom=chrom.replace(".", "_")) for chrom in regions.keys()
]
chunks_join = [
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
    .format(all_results=", ".join(map(repr, result_files_join)))
]

# Join results
tmpdir_join = tempfile.mkdtemp("snake_par")

resources = ResourceUsage(cores=1, memory=gib(1.0), duration=hours(72.0))
for i, (chrom, chrom_regions) in enumerate(regions.items()):
    chrom_name = chrom.replace(".", "_")
    chunks_join.append(
        textwrap.dedent(
            r"""
    rule chunk_{chrom}:
        output:
            touch('job_out.{chrom}.d/.done'),
        shell:
            r'''
            set -euo pipefail

            # Hack: get back bin directory of base/root environment.
            export PATH=$PATH:$(dirname $(dirname $(which conda)))/bin

            # concatenate files -------------------------------------------------------
            for jobno in {{{{0..{max_jobno}}}}}; do
                bgzip -@ 4 -c -d {work_dir}/job_out.{chrom}.$jobno.d/out/tmp_{chrom}.$jobno.g.vcf.gz
            done \
            | snappy-vcf_first_header \
            | bgzip -@ 4 -c \
            > job_out.{chrom}.d/{output_vcf}

            tabix -f job_out.{chrom}.d/{output_vcf}

            # compute MD5 sums
            pushd job_out.{chrom}.d && for f in *; do md5sum $f >$f.md5; done && popd

            # move out
            mkdir -p {output_dir}
            mv -f job_out.{chrom}.d/* {output_dir}
            '''

    cluster_config['chunk_{chrom}'] = {resources}
    """
        ).format(
            chrom=chrom_name,
            output_vcf=(
                "{snakemake.wildcards.mapper}.gatk_hc_gvcf." "combine_gvcf.{chrom}.g.vcf.gz"
            ).format(snakemake=snakemake, chrom=chrom_name),
            resources=repr(SgeResourceUsageConverter(resources).to_res_dict()),
            work_dir=tmpdir_work,
            output_dir=os.path.dirname(snakemake.output.vcf[0]),
            max_jobno=len(chrom_regions) - 1,
        )
    )

with in_working_dir(tmpdir_join, print_chdir=True):
    with open("Snakefile", "wt") as snakefile:
        print("\n\n".join(chunks_join), file=sys.stderr)
        print("\n\n".join(chunks_join), file=snakefile)

# Write out config file
with in_working_dir(tmpdir_join, print_chdir=True):
    with open("config.json", "wt") as configfile:
        json.dump(snakemake.config, configfile)
    with open("config.json", "rt") as configfile:
        print(configfile.read(), file=sys.stderr)

# Launch execution of Snakefile
with in_working_dir(tmpdir_join, print_chdir=True):
    run_snakemake(
        snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"],
        num_jobs=snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"].get(
            "num_jobs_combine_gvcf_cohort"
        ),
        job_name_token="gatk_hc_gvcf_par_combine_gvcf_join",
        drmaa_snippet=(
            snakemake.config["step_config"]["variant_calling"]["gatk_hc_gvcf"]["drmaa_snippet"]
            or snakemake.config["step_config"]["variant_calling"]["drmaa_snippet"]
        ),
    )

# -*- coding: utf-8 -*-
"""Definition for Mutect2 variant caller in parallel, genome is split into windows

isort:skip_file
"""

import os
import sys
import textwrap

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelVariantCallingBaseWrapper,
    gib_to_string,
    hours,
)


class ParallelMutect2Wrapper(ParallelVariantCallingBaseWrapper):
    """Parallel execution of MuTect 2"""

    inner_wrapper = "mutect2/prepare_panel"
    step_name = "panel_of_normals"
    tool_name = "mutect2"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=1,
            memory=gib_to_string(14.0 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=1,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )

    def construct_preamble(self):
        """Return a preamble that redefines resource_chunk_{threads,memory} to
        define functions as "scaling up" with the number of attempts.
        """
        return (
            textwrap.dedent(
                r"""
            shell.executable("/bin/bash")
            shell.prefix("set -ex;")

            configfile: 'config.json'

            localrules: all

            def multiply_time(day_time_str, factor):
                # Check if time contains day, ex: '1-00:00:00'
                if "-" in day_time_str:
                    arr_ = day_time_str.split("-")
                    days = int(arr_[0])
                    time_str = arr_[1]
                else:
                    days = 0
                    time_str = day_time_str

                # Process based on time structure
                arr_ = time_str.split(":")
                if time_str.count(":") == 2: # hours:minutes:seconds
                    seconds = int(arr_[0]) * 60 * 60 + int(arr_[1]) * 60 + int(arr_[2])
                elif time_str.count(":") == 1: # minutes:seconds
                    seconds = int(arr_[0]) * 60 + int(arr_[1])
                elif time_str.count(":") == 0: # minutes
                    seconds = int(time_str) * 60
                else:
                    raise ValueError(f"Invalid time: {{day_time_str}}")
                # Add days to second
                seconds += days * 86400

                # Apply factor
                seconds = int(seconds * factor)

                # Normalise time
                (norm_days, remainder) = divmod(seconds, 86400)
                (hours, remainder) = divmod(remainder, 3600)
                (minutes, seconds) = divmod(remainder, 60)

                # Fill string - example hour '7' -> '07'
                h_str = str(hours).zfill(2)
                m_str = str(minutes).zfill(2)
                s_str = str(seconds).zfill(2)

                return "%d-%s:%s:%s" % (norm_days, h_str, m_str, s_str)


            def multiply_memory(memory_str, factor):
                memory_mb = None
                suffixes = (
                    ("k", 1e-3),
                    ("M", 1),
                    ("G", 1e3),
                    ("T", 1e6),
                )
                for (suffix, mult) in suffixes:
                    if memory_str.endswith(suffix):
                        memory_mb = float(memory_str[:-1]) * mult
                        break
                # No match, assume no suffix int
                if not memory_mb:
                    memory_mb = float(memory_str)
                return int(memory_mb * factor)

            def resource_chunk_threads(wildcards):
                '''Return the number of threads to use for running one chunk.'''
                return {chunk_resources_threads}

            def resource_chunk_memory(wildcards, attempt):
                '''Return the memory to use for running one chunk.'''
                return multiply_memory({chunk_resources_memory}, attempt)

            def resource_chunk_time(wildcards, attempt):
                '''Return the time to use for running one chunk.'''
                return multiply_time({chunk_resources_time}, attempt)

            def resource_chunk_partition(wildcards):
                '''Return the partition to use for running one chunk.'''
                return {chunk_resources_partition}

            def resource_merge_threads(wildcards):
                '''Return the number of threads to use for running merging.'''
                return {merge_resources_threads}

            def resource_merge_memory(wildcards):
                '''Return the memory to use for running merging.'''
                return {merge_resources_memory}

            def resource_merge_time(wildcards):
                '''Return the time to use for running merging.'''
                return {merge_resources_time}

            def resource_merge_partition(wildcards):
                '''Return the partition to use for running merging.'''
                return {merge_resources_partition}

            rule all:
                input: **{all_output}
        """
            )
            .lstrip()
            .format(
                all_output=repr(self.get_all_output()),
                chunk_resources_threads=repr(self.job_resources.threads),
                chunk_resources_time=repr(self.job_resources.time),
                chunk_resources_memory=repr(self.job_resources.memory),
                chunk_resources_partition=repr(self.job_resources.partition),
                merge_resources_threads=repr(self.merge_resources.threads),
                merge_resources_time=repr(self.merge_resources.time),
                merge_resources_memory=repr(self.merge_resources.memory),
                merge_resources_partition=repr(self.merge_resources.partition),
            )
        )

    def _construct_level_one_merge_rule(self, chunk_no, merge_input):
        return (
            textwrap.dedent(
                r"""
            rule merge_chunk_{chunk_no}:
                input: {chunk_input}
                output:
                    vcf='merge_out.{chunk_no}.d/out/out.vcf.gz',
                    tbi='merge_out.{chunk_no}.d/out/out.vcf.gz.tbi',
                    stats='merge_out.{chunk_no}.d/out/out.vcf.stats',
                    f1r2='merge_out.{chunk_no}.d/out/out.f1r2_tar.tar.gz'
                threads: resource_chunk_threads
                resources:
                    time=resource_chunk_time,
                    memory=resource_chunk_memory,
                    partition=resource_chunk_partition,
                shell:
                    r'''
                    set -euo pipefail  # Unofficial Bash strict mode

                    # Concatenate VCF files -----------------------------------------------

                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o {{output.vcf}} \
                        -O z \
                        {{input}}

                    tabix -f {{output.vcf}}
                    '''
        """
            )
            .lstrip()
            .format(
                chunk_no=chunk_no,
                chunk_input=repr(merge_input),
            )
        )

    def _construct_final_merge_rule(self, merge_input):
        return (
            textwrap.dedent(
                r"""
            rule merge_all:
                input: {all_input}
                output: **{all_output}
                threads: resource_merge_threads
                resources:
                    time=resource_merge_time,
                    memory=resource_merge_memory,
                    partition=resource_merge_partition,
                log: **{all_log}
                shell:
                    r'''
                    set -euo pipefail  # Unofficial Bash strict mode

                   # Initialize output directory -----------------------------------------

                    outdir=$(basename {{output.vcf}})

                    mkdir -p output

                    # Concatenate VCF files -----------------------------------------------

                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o output/out.vcf.gz \
                        -O z \
                        {{input}}

                    # tabix index & md5 checksums -----------------------------------------

                    tabix -f output/out.vcf.gz

                    pushd output
                    for f in *; do md5sum $f >$f.md5; done
                    popd

                    # Move to output directory --------------------------------------------

                    mkdir -p $(dirname {{output.vcf}})
                    mv output/out.vcf.gz {{output.vcf}}
                    mv output/out.vcf.gz.md5 {{output.vcf_md5}}
                    mv output/out.vcf.gz.tbi {{output.tbi}}
                    mv output/out.vcf.gz.tbi.md5 {{output.tbi_md5}}

                    # Write out information about conda installation.
                    conda list >{{log.conda_list}}
                    conda info >{{log.conda_info}}

                    pushd $(dirname {{log.conda_list}})
                    md5sum $(basename {{log.conda_list}}) >$(basename {{log.conda_list}}).md5
                    md5sum $(basename {{log.conda_info}}) >$(basename {{log.conda_info}}).md5
                    popd
                    '''
        """
            )
            .lstrip()
            .format(
                all_input=repr(merge_input),
                all_output=repr(self.get_all_output()),
                all_log=repr(self.get_all_log_files()),
            )
        )

    def construct_parallel_rules(self):
        """Construct the rules for parallel processing to generate."""
        for jobno, region in enumerate(self.get_regions()):
            params = dict(self.snakemake.params)
            params.setdefault("args", {}).update({"intervals": [region.human_readable()]})
            output = {
                key: "job_out.{jobno}.d/out/tmp_{jobno}.{ext}".format(jobno=jobno, ext=ext)
                for key, ext in self.key_ext.items()
            }
            vals = {
                "input_bam": repr(self.snakemake.input.normal_bam),
                "jobno": jobno,
                "params": repr(params),
                "output": repr(output),
                "wrapper_prefix": "file://" + self.wrapper_base_dir,
                "inner_wrapper": self.inner_wrapper,
            }
            yield textwrap.dedent(
                r"""
                rule chunk_{jobno}:
                    input:
                        {input_bam},
                    output:
                        touch("job_out.{jobno}.d/.done"),
                        **{output}
                    threads: resource_chunk_threads
                    resources:
                        time=resource_chunk_time,
                        memory=resource_chunk_memory,
                        partition=resource_chunk_partition,
                    params:
                        **{params}
                    wrapper: '{wrapper_prefix}/snappy_wrappers/wrappers/{inner_wrapper}'
            """
            ).format(**vals).lstrip()

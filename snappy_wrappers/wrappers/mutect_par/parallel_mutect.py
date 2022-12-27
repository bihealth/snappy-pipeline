# -*- coding: utf-8 -*-
"""Definition for Mutect variant caller in parallel, genome is split into windows
"""

import os
import sys
import textwrap

# The following is required for being able to import snappy_wrappers modules
# inside wrappers.  These run in an "inner" snakemake process which uses its
# own conda environment which cannot see the snappy_pipeline installation.
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.resource_usage import ResourceUsage  # noqa: E402
from snappy_wrappers.wrapper_parallel import (  # noqa: E402
    ParallelSomaticVariantCallingBaseWrapper,
    gib_to_string,
    hours,
)


class ParallelMutectWrapper(ParallelSomaticVariantCallingBaseWrapper):
    """Parallel execution of MuTect"""

    # TODO: probably, nobody looked at anything but the vcf/tbi files... get rid of them?
    realpath_output_keys = (
        "vcf",
        "vcf_md5",
        "vcf_tbi",
        "vcf_tbi_md5",
        "full_vcf",
        "full_vcf_md5",
        "full_vcf_tbi",
        "full_vcf_tbi",
        "full_vcf_tbi_md5",
        "txt",
        "txt_md5",
        "wig",
        "wig_md5",
    )

    key_ext = {
        "txt": "full.out.txt.gz",
        "txt_md5": "full.out.txt.gz.md5",
        "vcf": "vcf.gz",
        "vcf_md5": "vcf.gz.md5",
        "vcf_tbi": "vcf.gz.tbi",
        "vcf_tbi_md5": "vcf.gz.tbi.md5",
        "full_vcf": "full.vcf.gz",
        "full_vcf_md5": "full.vcf.gz.md5",
        "full_vcf_tbi": "full.vcf.gz.tbi",
        "full_vcf_tbi_md5": "full.vcf.gz.tbi.md5",
        "wig": "full.wig.txt.gz",
        "wig_md5": "full.wig.txt.gz.md5",
    }

    inner_wrapper = "mutect"
    step_name = "somatic_variant_calling"
    tool_name = "mutect"

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            threads=1,
            memory=gib_to_string(7.5 * self.get_job_mult_memory()),
            time=hours(4 * self.get_job_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )
        self.merge_resources = ResourceUsage(
            threads=1,
            memory=gib_to_string(2.0 * self.get_merge_mult_memory()),
            time=hours(4 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )

    def construct_merge_rule(self):
        """Join the overall result files"""
        return (
            textwrap.dedent(
                r"""
            rule merge_all:
                input: [{all_input}]
                output: **{all_output}
                threads: resource_merge_threads
                resources:
                    time=resource_merge_time,
                    memory=resource_merge_memory,
                    partition=resource_merge_partition,
                log: **{all_log}
                shell:
                    r'''
                    # Initialize output directory -----------------------------------------

                    outdir=$(basename {{output.vcf}})

                    mkdir -p output

                    # Take first header -------------------------------------------------------

                    set +o pipefail
                    zcat job_out.0.d/out/tmp_0.full.out.txt.gz | head -n 2 > output/result.full.out.txt
                    zcat job_out.0.d/out/tmp_0.full.wig.txt.gz | head -n 2 > output/result.full.wig.txt
                    set -o pipefail

                    # Append body contents ----------------------------------------------------

                    for jobno in {{{{0..{max_jobno}}}}}; do
                        set +o pipefail
                        zcat job_out.$jobno.d/out/tmp_$jobno.full.out.txt.gz | tail -n +3 >> output/result.full.out.txt
                        zcat job_out.$jobno.d/out/tmp_$jobno.full.wig.txt.gz | tail -n +3 >> output/result.full.wig.txt
                        set -o pipefail
                    done

                    # Use bcftools concat for VCF files ---------------------------------------

                    bcftools concat -a -d none -O z -o output/result.full.vcf.gz job_out.*.d/out/tmp_*.full.vcf.gz
                    bcftools concat -a -d none -O z -o output/result.vcf.gz job_out.*.d/out/tmp_*[0-9].vcf.gz

                    # bgzip output and create tabix index -------------------------------------

                    bgzip -f output/result.full.out.txt
                    bgzip -f output/result.full.wig.txt

                    tabix -f output/result.full.vcf.gz
                    tabix -f output/result.vcf.gz

                    pushd output
                    for f in *; do md5sum $f >$f.md5; done
                    popd

                    # Move to output directory ------------------------------------------------

                    mkdir -p $(dirname {{output.txt}})
                    mv output/result.full.out.txt.gz {{output.txt}}
                    mv output/result.full.out.txt.gz.md5 {{output.txt_md5}}
                    mv output/result.full.vcf.gz {{output.full_vcf}}
                    mv output/result.full.vcf.gz.md5 {{output.full_vcf_md5}}
                    mv output/result.full.vcf.gz.tbi {{output.full_vcf_tbi}}
                    mv output/result.full.vcf.gz.tbi.md5 {{output.full_vcf_tbi_md5}}
                    mv output/result.vcf.gz {{output.vcf}}
                    mv output/result.vcf.gz.md5 {{output.vcf_md5}}
                    mv output/result.vcf.gz.tbi {{output.vcf_tbi}}
                    mv output/result.vcf.gz.tbi.md5 {{output.vcf_tbi_md5}}
                    mv output/result.full.wig.txt.gz {{output.wig}}
                    mv output/result.full.wig.txt.gz.md5 {{output.wig_md5}}

                    # Write out information about conda installation.
                    conda list >{{log.conda_list}}
                    conda info >{{log.conda_info}}
                    md5sum {{log.conda_list}} >{{log.conda_list_md5}}
                    md5sum {{log.conda_info}} >{{log.conda_info_md5}}
                    '''
        """
            )
            .lstrip()
            .format(
                all_input=", ".join(map(repr, self.construct_parallel_result_files())),
                all_output=repr(self.get_all_output()),
                all_log=repr(self.get_all_log_files()),
                max_jobno=len(self.get_regions()) - 1,
            )
        )

# -*- coding: utf-8 -*-
"""Wrapper for running MuTect 2 variant caller in parallel, genome is split into windows
"""

import os
import sys
import textwrap

from snakemake import shell

# A hack is required for being able to import snappy_wrappers modules when in development mode.
# TODO: is there a more elegant way?
base_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, base_dir)

from snappy_wrappers.wrapper_parallel import (
    ParallelSomaticVariantCallingBaseWrapper,
    ResourceUsage,
    gib,
    hours,
)


class ParallelMutect2Wrapper(ParallelSomaticVariantCallingBaseWrapper):
    """Parallel execution of MuTect 2"""

    inner_wrapper = "mutect2/run"
    step_name = "somatic_variant_calling"
    tool_name = "mutect2"

    realpath_output_keys = (
        "raw",
        "raw_md5",
        "raw_tbi",
        "raw_tbi_md5",
        "stats",
        "stats_md5",
        "f1r2",
        "f1r2_md5",
    )
    key_ext = {
        "raw": "vcf.gz",
        "raw_md5": "vcf.gz.md5",
        "raw_tbi": "vcf.gz.tbi",
        "raw_tbi_md5": "vcf.gz.tbi.md5",
        "stats": "vcf.stats",
        "stats_md5": "vcf.stats.md5",
        "f1r2": "f1r2_tar.tar.gz",
        "f1r2_md5": "f1r2_tar.tar.gz.md5",
    }

    def __init__(self, snakemake):
        super().__init__(snakemake)
        self.job_resources = ResourceUsage(
            cores=1,
            memory=gib(14.0 * self.get_job_mult_memory()),
            duration=hours(1 * self.get_job_mult_time()),
        )
        self.merge_resources = ResourceUsage(
            cores=1,
            memory=gib(16.0 * self.get_merge_mult_memory()),
            duration=hours(1 * self.get_merge_mult_time()),
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
                shell:
                    r'''
                    set -euo pipefail  # inofficial Bash strict mode

                    # Concatenate VCF files -----------------------------------------------

                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o {{output.vcf}} \
                        -O z \
                        {{input}}

                    tabix -f {{output.vcf}}

                    # Concatenate stats ---------------------------------------------------

                    chunks=$(echo "{{input}}" | sed -e "s/\.vcf\.gz/.vcf.stats/g" | sed -e "s/ / -stats /")
                    gatk MergeMutectStats -stats $chunks -O {output.stats}

                    # Concatenate f1r2 tar files ------------------------------------------

                    tar_dir=$(dirname "{{output.vcf}}")
                    tar_dir="${{tar_dir}}/out.f1r2_tar"

                    mkdir -p $tar_dir
                    pushd $tar_dir

                    chunks=$(echo "{{input}}" | sed -e "s/\.vcf\.gz/.f1r2_tar.tar.gz/" | sed -re "s/job_out\.([0-9]+)\.d/..\/..\/job_out.\1.d/g")
                    for chunk in $chunks
                    do
                        tar -zxvf $chunk
                    done
                    tar -zcvf ../out.f1r2_tar.tar.gz *

                    popd
                    rm -rf $tar_dir
                    '''

            cluster_config['merge_chunk_{chunk_no}'] = {resources}
        """
            )
            .lstrip()
            .format(
                chunk_no=chunk_no,
                chunk_input=repr(merge_input),
                resources=repr(self.res_converter(self.merge_resources).to_res_dict()),
            )
        )

    def _construct_final_merge_rule(self, merge_input):
        return (
            textwrap.dedent(
                r"""
            rule merge_all:
                input: {all_input}
                output: **{all_output}
                log: **{all_log}
                shell:
                    r'''
                    set -euo pipefail  # inofficial Bash strict mode

                    # Initialize output directory -----------------------------------------

                    outdir=$(basename {{output.raw}})

                    mkdir -p output

                    # Concatenate VCF files -----------------------------------------------

                    bcftools concat \
                        --allow-overlaps \
                        -d none \
                        -o output/out.vcf.gz \
                        -O z \
                        {{input}}

                    # Concatenate stats ---------------------------------------------------

                    chunks=$(echo "{{input}}" | sed -e "s/\.vcf\.gz/.vcf.stats/g" | sed -e "s/ / -stats /g")
                    stats="output/out.vcf.stats"
                    gatk MergeMutectStats -stats $chunks -O $stats

                    # Concatenate f1r2 tar files ------------------------------------------

                    tar_dir="output/out.f1r2_tar"

                    mkdir -p $tar_dir
                    pushd $tar_dir

                    chunks=$(echo "{{input}}" | sed -e "s/\.vcf\.gz/.f1r2_tar.tar.gz/g" | sed -re "s/(job|merge)_out\.([0-9]+)\.d/..\/..\/\1_out.\2.d/g")
                    for chunk in $chunks
                    do
                        tar -zxvf $chunk
                    done
                    tar -zcvf ../out.f1r2_tar.tar.gz *

                    popd
                    rm -rf $tar_dir

                    # tabix index & md5 checksums -----------------------------------------

                    tabix -f output/out.vcf.gz

                    pushd output
                    for f in *; do md5sum $f >$f.md5; done
                    popd

                    # Move to output directory --------------------------------------------

                    mkdir -p $(dirname {{output.raw}})
                    mv output/out.vcf.gz {{output.raw}}
                    mv output/out.vcf.gz.md5 {{output.raw_md5}}
                    mv output/out.vcf.gz.tbi {{output.raw_tbi}}
                    mv output/out.vcf.gz.tbi.md5 {{output.raw_tbi_md5}}
                    mv output/out.vcf.stats {{output.stats}}
                    mv output/out.vcf.stats.md5 {{output.stats_md5}}
                    mv output/out.f1r2_tar.tar.gz {{output.f1r2}}
                    mv output/out.f1r2_tar.tar.gz.md5 {{output.f1r2_md5}}

                    # Write out information about conda installation.
                    conda list >{{log.conda_list}}
                    conda info >{{log.conda_info}}

                    pushd $(dirname {{log.conda_list}})
                    md5sum $(basename {{log.conda_list}}) >$(basename {{log.conda_list}}).md5
                    md5sum $(basename {{log.conda_info}}) >$(basename {{log.conda_info}}).md5
                    popd
                    '''

            cluster_config['merge_all'] = {resources}
        """
            )
            .lstrip()
            .format(
                all_input=repr(merge_input),
                all_output=repr(self.get_all_output()),
                all_log=repr(self.get_all_log_files()),
                resources=repr(self.res_converter(self.merge_resources).to_res_dict()),
            )
        )


# Kick off execution using the wrapper class defined above.
ParallelMutect2Wrapper(snakemake).run().shutdown_logging()

# Compute MD5 sums of logs.
shell(
    r"""
md5sum {snakemake.log.log} >{snakemake.log.log_md5}
"""
)

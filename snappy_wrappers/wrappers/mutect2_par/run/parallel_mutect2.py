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
    gib_to_string,
    hours,
    ParallelMutect2BaseWrapper,
)


class ParallelMutect2Wrapper(ParallelMutect2BaseWrapper):

    inner_wrapper = "mutect2/run"
    step_name = "somatic_variant_calling"
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
            memory=gib_to_string(16.0 * self.get_merge_mult_memory()),
            time=hours(1 * self.get_merge_mult_time()),
            partition=os.getenv("SNAPPY_PIPELINE_PARTITION"),
        )

    @classmethod
    def allow_resources_increase(cls):
        return True

    @classmethod
    def merge_code_level_one(cls):
        return cls._merge_code()

    @classmethod
    def merge_code_final(cls):
        return cls._merge_code()

    @classmethod
    def _merge_code(cls):
        return textwrap.dedent(
            r"""
        # Concatenate raw calls vcfs & index result ----------------------
        bcftools concat \
            --allow-overlaps \
            -d none \
            -o {output.raw} \
            -O z \
            {input.raw}
        tabix -f {output.raw}

        # Concatenate stats with GATK tool -------------------------------
        stats=$(echo "{input.stats}" | sed -e "s/ / -stats /g")
        gatk MergeMutectStats -stats $stats -O {output.stats}

        # Contatenate orientation tar files ------------------------------
        tmpdir=$(mktemp -d)
        for tar_file in {input.f1r2}
        do
            abs_path=$(realpath $tar_file)
            pushd $tmpdir
            tar -zxvf $abs_path
            popd
        done
        tar -zcvf {output.f1r2} -C $tmpdir .
        rm -rf $tmpdir

        # Compute md5 sums -----------------------------------------------
        pushd $(dirname {output.raw})
        for f in *; do md5sum $f > $f.md5; done
        popd
                """
        )

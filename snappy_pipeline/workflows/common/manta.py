"""Workflow step parts for Manta.

These are used in both ``sv_calling_targeted`` and ``sv_calling_wgs``.
"""

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart
from snappy_pipeline.workflows.common.sv_calling import (
    VCF_OUTPUTS,
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_wrappers.resource_usage import ResourceUsage


class MantaStepPart(
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Manta"""

    name = "manta"
    actions = ("run",)

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config.manta.num_threads
        return ResourceUsage(
            threads=num_threads,
            runtime="7d",  # 7 days
            mem=f"{int(3.5 * 1024 * num_threads)}MB",
        )

    @dictify
    def _get_input_files_run(self, wildcards):
        bams = [
            self.parent.get_upstream_paths("alignments", library_name=library_name).bam
            for library_name in self.parent.pedigrees(include_background=True)[
                wildcards.library_name
            ]
        ]
        yield "bam", bams
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    def _get_output_files_run(self):
        return dict(VCF_OUTPUTS)

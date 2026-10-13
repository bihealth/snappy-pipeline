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

    def __init__(self, parent):
        super().__init__(parent)
        #: Shortcuts from index NGS library name to Pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_run(self, wildcards):
        bams = []
        for donor in self.index_ngs_library_to_pedigree[wildcards.library_name].donors:
            if donor.dna_ngs_library:
                alignments = self.parent.get_upstream_paths(
                    "alignments", library_name=donor.dna_ngs_library.name
                )
                bams.append(alignments.bam)
        yield "bam", bams
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    def _get_output_files_run(self):
        return dict(VCF_OUTPUTS)

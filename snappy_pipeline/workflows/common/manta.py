"""Workflow step parts for Manta.

These are used in both ``sv_calling_targeted`` and ``sv_calling_wgs``.
"""

from typing import Any

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart
from snappy_pipeline.workflows.abstract.common import (
    ForwardResourceUsageMixin,
    ForwardSnakemakeFilesMixin,
    augment_work_dir_with_output_links,
)
from snappy_pipeline.workflows.common.sv_calling import (
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_wrappers.resource_usage import ResourceUsage


class MantaStepPart(
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Manta"""

    name = "manta"
    actions = ("run",)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        if action not in self.actions:
            actions_str = ", ".join(self.actions)
            error_message = f"Action '{action}' is not supported. Valid options: {actions_str}"
            raise UnsupportedActionException(error_message)
        num_threads = self.config.manta.num_threads
        return ResourceUsage(
            threads=num_threads,
            time="7-00:00:00",  # 3 days
            memory=f"{int(3.5 * 1024 * num_threads)}M",
        )

    def __init__(self, parent):
        super().__init__(parent)
        #: Shortcuts from index NGS library name to Pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_run(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bams = []
        for donor in self.index_ngs_library_to_pedigree[wildcards.library_name].donors:
            if donor.dna_ngs_library:
                token = f"{wildcards.mapper}.{donor.dna_ngs_library.name}"
                bams.append(ngs_mapping(f"output/{token}/out/{token}.bam"))
        yield "bam", bams

    @dictify
    def _get_output_files_run(self):
        infix = "{mapper}.manta.{library_name}"
        work_files = {
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from augment_work_dir_with_output_links(
            work_files, self.get_log_file().values()
        ).items()
    
    def get_args(self, action: str) -> dict[str, Any]:
        self._validate_action(action)
        return {"reference": self.parent.w_config.static_data_config.reference.path}
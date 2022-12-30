"""Workflow step parts for Manta.

These are used in both ``sv_calling_targeted`` and ``sv_calling_wgs``.
"""

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
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Manta"""

    name = "manta"
    actions = ("run",)

    resource_usage_dict = {
        "run": ResourceUsage(
            threads=16,
            time="2-00:00:00",
            memory=f"{int(3.75 * 1024 * 16)}M",
        )
    }

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
        infix = "{mapper}.delly2.{library_name}"
        work_files = {
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from augment_work_dir_with_output_links(
            work_files, self.get_log_file("merge_genotypes").values()
        ).items()

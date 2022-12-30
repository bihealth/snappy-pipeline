"""Workflow step parts for Delly.

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


class Delly2StepPart(
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Delly2"""

    name = "delly2"
    actions = ("call", "merge_calls", "genotype", "merge_genotypes")

    _cheap_resource_usage = ResourceUsage(
        threads=2,
        time="4-00:00:00",
        memory=f"{7 * 1024 * 2}M",
    )
    _normal_resource_usage = ResourceUsage(
        threads=2,
        time="7-00:00:00",  # 7 days
        memory=f"{20 * 1024 * 2}M",
    )
    resource_usage_dict = {
        "call": _normal_resource_usage,
        "merge_calls": _cheap_resource_usage,
        "genotype": _normal_resource_usage,
        "merge_genotypes": _cheap_resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)

        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

        self.donor_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_call(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        token = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{token}/out/{token}.bam")

    @dictify
    def _get_output_files_call(self):
        infix = "{mapper}.delly2_call.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_calls(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.delly2_call.{donor.dna_ngs_library.name}"
                bcfs.append(f"work/{infix}/out/{infix}.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_calls(self):
        infix = "{mapper}.delly2_merge_calls.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield from self._get_input_files_call(wildcards).items()
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        infix = f"{wildcards.mapper}.delly2_merge_calls.{pedigree.index.dna_ngs_library.name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"

    @dictify
    def _get_output_files_genotype(self):
        infix = "{mapper}.delly2_genotype.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_genotypes(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.delly2_genotype.{donor.dna_ngs_library.name}"
                bcfs.append(f"work/{infix}/out/{infix}.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_genotypes(self):
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

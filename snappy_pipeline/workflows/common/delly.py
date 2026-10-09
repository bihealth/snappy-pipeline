"""Workflow step parts for Delly.

These are used in both ``sv_calling_targeted`` and ``sv_calling_wgs``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

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
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Delly2"""

    name = "delly2"
    actions = ("call", "merge_calls", "genotype", "merge_genotypes")

    _cheap_resource_usage = ResourceUsage(
        threads=1,
        runtime="1d",
        mem="4GB",
    )
    _normal_resource_usage = ResourceUsage(threads=1, runtime="2d", mem="16GB")
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

    def _get_params_call(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "genome": self.parent.get_upstream_paths("reference").fasta,
            "config": dict(self.config.get(self.name)),
        }

    _get_params_merge_calls = _get_params_call
    _get_params_genotype = _get_params_call
    _get_params_merge_genotypes = _get_params_call

    @dictify
    def _get_input_files_call(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai

    @dictify
    def _get_output_files_call(self):
        prefix = "work/{library_name}/out/{library_name}.call"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_md5", f"{prefix}.bcf.md5"
        yield "bcf_csi", f"{prefix}.bcf.csi"
        yield "bcf_csi_md5", f"{prefix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_calls(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                library_name = donor.dna_ngs_library.name
                bcfs.append(f"work/{library_name}/out/{library_name}.call.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_calls(self):
        prefix = "work/{library_name}/out/{library_name}.merge_calls"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_md5", f"{prefix}.bcf.md5"
        yield "bcf_csi", f"{prefix}.bcf.csi"
        yield "bcf_csi_md5", f"{prefix}.bcf.csi.md5"

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield from self._get_input_files_call(wildcards).items()
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        index_library_name = pedigree.index.dna_ngs_library.name
        yield "bcf", f"work/{index_library_name}/out/{index_library_name}.merge_calls.bcf"

    @dictify
    def _get_output_files_genotype(self):
        prefix = "work/{library_name}/out/{library_name}.genotype"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_md5", f"{prefix}.bcf.md5"
        yield "bcf_csi", f"{prefix}.bcf.csi"
        yield "bcf_csi_md5", f"{prefix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_genotypes(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                library_name = donor.dna_ngs_library.name
                bcfs.append(f"work/{library_name}/out/{library_name}.genotype.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_genotypes(self):
        prefix = "work/{library_name}/out/{library_name}"
        work_files = {
            "vcf": f"{prefix}.vcf.gz",
            "vcf_md5": f"{prefix}.vcf.gz.md5",
            "vcf_tbi": f"{prefix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"{prefix}.vcf.gz.tbi.md5",
        }
        yield from augment_work_dir_with_output_links(
            work_files, self.get_log_file("merge_genotypes").values()
        ).items()

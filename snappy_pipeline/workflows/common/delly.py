"""Workflow step parts for Delly.

These are used in both ``sv_calling_targeted`` and ``sv_calling_wgs``.
"""

from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart
from snappy_pipeline.workflows.common.sv_calling import (
    VCF_OUTPUTS,
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_wrappers.resource_usage import ResourceUsage


class Delly2StepPart(
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Delly2"""

    name = "delly2"
    actions = ("call", "merge_calls", "genotype", "merge_genotypes")

    #: ``genotype`` runs ``delly call`` again, but takes the arguments of the discovery only
    extra_arguments_action = "call"

    intermediate_outputs = {
        "call": ("bcf", "bcf_csi"),
        "merge_calls": ("bcf", "bcf_csi"),
        "genotype": ("bcf", "bcf_csi"),
    }

    _cheap_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="4GB")
    _normal_resource_usage = ResourceUsage(threads=1, runtime="2d", mem="16GB")
    resource_usage = {
        "call": _normal_resource_usage,
        "merge_calls": _cheap_resource_usage,
        "genotype": _normal_resource_usage,
        "merge_genotypes": _cheap_resource_usage,
    }

    def _get_params_call(self, wildcards: Wildcards) -> dict[str, Any]:
        config = dict(self.config.get(self.name))
        return {
            "config": {
                key: config[key] for key in ("map_qual", "qual_tra", "geno_qual", "mad_cutoff")
            },
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
        yield from self.parent.reference_inputs().items()
        if self.config.depends_on.excluded_regions:
            yield "exclude", self.parent.get_upstream_paths("excluded_regions").bed

    @dictify
    def _get_output_files_call(self):
        prefix = "work/{library_name}/out/{library_name}.call"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_csi", f"{prefix}.bcf.csi"

    @dictify
    def _get_input_files_merge_calls(self, wildcards):
        bcfs = [
            f"work/{library_name}/out/{library_name}.call.bcf"
            for library_name in self.parent.pedigrees(include_background=True)[
                wildcards.library_name
            ]
        ]
        yield "bcf", bcfs
        yield "bcf_csi", [f"{bcf}.csi" for bcf in bcfs]

    @dictify
    def _get_output_files_merge_calls(self):
        prefix = "work/{library_name}/out/{library_name}.merge_calls"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_csi", f"{prefix}.bcf.csi"

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield from self._get_input_files_call(wildcards).items()
        index_library_name = self.parent.pedigree_index(
            wildcards.library_name, include_background=True
        )
        bcf = f"work/{index_library_name}/out/{index_library_name}.merge_calls.bcf"
        yield "bcf", bcf
        yield "bcf_csi", f"{bcf}.csi"

    @dictify
    def _get_output_files_genotype(self):
        prefix = "work/{library_name}/out/{library_name}.genotype"
        yield "bcf", f"{prefix}.bcf"
        yield "bcf_csi", f"{prefix}.bcf.csi"

    @dictify
    def _get_input_files_merge_genotypes(self, wildcards):
        bcfs = [
            f"work/{library_name}/out/{library_name}.genotype.bcf"
            for library_name in self.parent.pedigrees(include_background=True)[
                wildcards.library_name
            ]
        ]
        yield "bcf", bcfs
        yield "bcf_csi", [f"{bcf}.csi" for bcf in bcfs]

    def _get_output_files_merge_genotypes(self):
        return dict(VCF_OUTPUTS)

# -*- coding: utf-8 -*-
"""Implementation of the gCNV common methods."""

import pandas as pd

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage


#: The kit name of whole genome models
WGS_KIT = "wgs"


def build_ngs_library_to_kit(
    libraries: pd.DataFrame, intervals: str, library_kits: dict[str, str]
) -> dict[str, str]:
    """Return the gCNV kit of each library in ``libraries``, the primary DNA library rows.

    ``intervals`` is the ``intervals`` field of the ``gcnv_model`` task.  With ``wgs`` every
    library gets the kit ``wgs``.  With ``targeted`` a library gets its sample sheet kit from
    ``library_kits`` (library name to kit), and libraries without a kit (not targeted libraries)
    are left out.
    """
    if intervals == "wgs":
        return {name: WGS_KIT for name in libraries["library_name"]}
    return {name: library_kits[name] for name in libraries["library_name"] if name in library_kits}


class GcnvWarning(UserWarning):
    """Base class for custom warning of the gCNV integration"""


class InconsistentLibraryKitsWarning(GcnvWarning):
    """Raised when library kits are not consistent within a pedigree

    Models can only be trained reliably on similar data.  This implies the same library kit.
    In the case that inconsistent library kits are used within a pedigree, this warning is
    raised.
    """


class TooFewSamplesWarning(GcnvWarning):
    """Raised when too few samples were provided"""


class PreprocessIntervalsCommonMixin:
    """Mixin used for the ``preprocess_intervals`` step."""

    def _get_input_files_preprocess_intervals(self, wildcards):
        inputs = self.parent.reference_inputs()
        # Targeted sequencing: the target regions of the library kit
        if self.config.intervals == "targeted":
            inputs["target_bed"] = self.parent.get_upstream_paths(
                "kits", library_kit=wildcards.library_kit
            ).targets
        return inputs

    @dictify
    def _get_output_files_preprocess_intervals(self):
        yield "interval_list", "work/{library_kit}/out/{library_kit}.interval_list"

    def _get_log_file_preprocess_intervals(self):
        return {"log": "work/{library_kit}/log/{library_kit}.preprocess_intervals.log"}


class CoverageCommonMixin:
    """Mixin used for ``coverage`` step"""

    @dictify
    def _get_input_files_coverage(self, wildcards):
        """Yield input files for ``coverage`` rule

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_name' (e.g., 'P001-N1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards
        """
        # Yield .interval list file.
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        yield "interval_list", self._interval_list_path(library_kit)
        # Yield input BAM and BAI files
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_output_files_coverage(self):
        yield "tsv", "work/{library_name}/out/{library_name}.coverage.tsv"

    def _get_log_file_coverage(self):
        return {"log": "work/{library_name}/log/{library_name}.coverage.log"}


class ContigPloidyCommonMixin:
    """Mixin used for ``contig_ploidy`` step"""

    def _get_log_file_contig_ploidy(self):
        return {"log": "work/{library_kit}/log/{library_kit}.contig_ploidy.log"}


class CallCnvsCommonMixin:
    """Mixin used for the ``call_cnvs`` step"""

    def _get_log_file_call_cnvs(self):
        return {"log": "work/{library_kit}/log/{library_kit}.{shard}.call_cnvs.log"}


class GcnvPostGermlineCallsCommonMixin:
    """Mixin used for the ``gcnv_post_germline_calls`` step"""

    def _get_log_file_post_germline_calls(self):
        return {"log": "work/{library_name}/log/{library_name}.post_germline_calls.log"}


class GcnvCommonStepPart(
    PreprocessIntervalsCommonMixin,
    CoverageCommonMixin,
    ContigPloidyCommonMixin,
    CallCnvsCommonMixin,
    GcnvPostGermlineCallsCommonMixin,
    BaseStepPart,
):
    """Class contains methods that are common for both gCNV model build and run."""

    #: Step name
    name = "gcnv"

    #: Dictionary: Key: library name (str); Value: library kit (str).
    ngs_library_to_kit = None

    def _interval_list_path(self, library_kit: str) -> str:
        """Return the preprocessed interval list of the kit, a local path of the model task."""
        return f"work/{library_kit}/out/{library_kit}.interval_list"

    #: Resource usage of the actions that run the gCNV model
    resource_usage = dict.fromkeys(
        ("call_cnvs", "post_germline_calls", "joint_germline_cnv_segmentation"),
        ResourceUsage(threads=16, runtime="4d", mem="46080MB"),
    )

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="7680MB")

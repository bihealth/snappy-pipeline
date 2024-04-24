# -*- coding: utf-8 -*-
"""Implementation of the gCNV common methods."""

from collections import OrderedDict

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage


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
        _ = wildcards
        return {}

    @dictify
    def _get_output_files_preprocess_intervals(self):
        ext = "interval_list"
        name_pattern = "gcnv_preprocess_intervals.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    def _get_log_file_preprocess_intervals(self):
        name_pattern = "gcnv_preprocess_intervals.{library_kit}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


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
        ext = "interval_list"
        library_kit = self.ngs_library_to_kit[wildcards.library_name]
        name_pattern = f"gcnv_preprocess_intervals.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"
        # Yield input BAM and BAI files
        ngs_mapping = self.parent.modules["ngs_mapping"]
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    @dictify
    def _get_output_files_coverage(self):
        ext = "tsv"
        name_pattern = "{mapper}.gcnv_coverage.{library_name}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    def _get_log_file_coverage(self):
        name_pattern = "{mapper}.gcnv_coverage.{library_name}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class ContigPloidyCommonMixin:
    """Mixin used for ``contig_ploidy`` step"""

    def _get_log_file_contig_ploidy(self):
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class CallCnvsCommonMixin:
    """Mixin used for the ``call_cnvs`` step"""

    def _get_log_file_call_cnvs(self):
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class GcnvPostGermlineCallsCommonMixin:
    """Mixin used for the ``gcnv_post_germline_calls`` step"""

    def _get_log_file_post_germline_calls(self):
        name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


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

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "high_resource": ResourceUsage(
            threads=16,
            time="4-00:00:00",
            memory="46080M",
        ),
        "default": ResourceUsage(
            threads=1,
            time="1-00:00:00",
            memory="7680M",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to donor
        self.index_ngs_library_to_donor = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_donor.update(sheet.index_ngs_library_to_donor)
        # Build shortcut from index library name to pedigree
        self.donor_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action):
        """Return input function for gCNV rule

        :param action: Action (i.e., step) in the workflow, examples: 'filter_intervals',
        'coverage'.
        :type action: str

        :return: Returns input function for gCNV rule based on inputted action.
        """
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def get_output_files(self, action):
        """Get output function for gCNV build model rule.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns output function for gCNV rule based on inputted action.
        """
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def get_log_file(self, action):
        """Get log file.

        :param action: Action (i.e., step) in the workflow, examples: 'filter_intervals',
        'coverage'.
        :type action: str

        :return: Returns template path to log file.
        """
        # TODO: move out the _get_log_file* functions where they belong
        self._validate_action(action)
        return getattr(self, f"_get_log_file_{action}")()

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action)
        high_resource_action_list = (
            "call_cnvs",
            "post_germline_calls",
            "joint_germline_cnv_segmentation",
        )
        if action in high_resource_action_list:
            return self.resource_usage_dict.get("high_resource")
        else:
            return self.resource_usage_dict.get("default")

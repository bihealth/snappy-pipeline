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
        inputs = {"reference": self.parent.get_upstream_paths("reference").fasta}
        # Targeted sequencing: the target regions of the library kit
        for item in getattr(self.config.gcnv, "path_target_interval_list_mapping", None) or []:
            if item.name == wildcards.library_kit:
                inputs["target_bed"] = item.path
                break
        return inputs

    @dictify
    def _get_output_files_preprocess_intervals(self):
        yield "interval_list", "work/{library_kit}/out/{library_kit}.interval_list"

    def _get_log_file_preprocess_intervals(self):
        return "work/{library_kit}/log/{library_kit}.preprocess_intervals.log"


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
        yield "interval_list", f"work/{library_kit}/out/{library_kit}.interval_list"
        # Yield input BAM and BAI files
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    @dictify
    def _get_output_files_coverage(self):
        yield "tsv", "work/{library_name}/out/{library_name}.coverage.tsv"

    def _get_log_file_coverage(self):
        return "work/{library_name}/log/{library_name}.coverage.log"


class ContigPloidyCommonMixin:
    """Mixin used for ``contig_ploidy`` step"""

    def _get_log_file_contig_ploidy(self):
        return "work/{library_kit}/log/{library_kit}.contig_ploidy.log"


class CallCnvsCommonMixin:
    """Mixin used for the ``call_cnvs`` step"""

    def _get_log_file_call_cnvs(self):
        return "work/{library_kit}/log/{library_kit}.{shard}.call_cnvs.log"


class GcnvPostGermlineCallsCommonMixin:
    """Mixin used for the ``gcnv_post_germline_calls`` step"""

    def _get_log_file_post_germline_calls(self):
        return "work/{library_name}/log/{library_name}.post_germline_calls.log"


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
            runtime="4d",
            mem="46080MB",
        ),
        "default": ResourceUsage(
            threads=1,
            runtime="1d",
            mem="7680MB",
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

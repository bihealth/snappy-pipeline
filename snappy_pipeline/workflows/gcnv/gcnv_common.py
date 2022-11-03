# -*- coding: utf-8 -*-
"""Implementation of the gCNV common methods.
"""

from collections import OrderedDict

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage


class GcnvStepPart(BaseStepPart):
    """Class contains methods that are common for both gCNV model build and run."""

    #: Step name
    name = "gcnv"

    #: Dictionary: Key: library name (str); Value: library kit (str).
    ngs_library_to_kit = None

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "high_resource": ResourceUsage(
            threads=16,
            time="2-00:00:00",
            memory="46080M",
        ),
        "default": ResourceUsage(
            threads=1,
            time="04:00:00",
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
        'coverage', 'extract_ped'.
        :type action: str

        :return: Returns input function for gCNV rule based on inputted action.
        """
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @staticmethod
    def _get_input_files_preprocess_intervals(wildcards):
        _ = wildcards
        return {}

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
        name_pattern = "gcnv_preprocess_intervals.{library_kit}".format(library_kit=library_kit)
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )
        # Yield input BAM and BAI files
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}{ext}"
        for key, ext in {"bam": ".bam", "bai": ".bam.bai"}.items():
            yield key, ngs_mapping(bam_tpl.format(ext=ext, **wildcards))

    def get_output_files(self, action):
        """Get output function for gCNV build model rule.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns output function for gCNV rule based on inputted action.
        """
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @staticmethod
    @dictify
    def _get_output_files_preprocess_intervals():
        ext = "interval_list"
        name_pattern = "gcnv_preprocess_intervals.{library_kit}"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @staticmethod
    @dictify
    def _get_output_files_coverage():
        ext = "tsv"
        name_pattern = "{mapper}.gcnv_coverage.{library_name}"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    def get_log_file(self, action):
        """Get log file.

        :param action: Action (i.e., step) in the workflow, examples: 'filter_intervals',
        'coverage', 'extract_ped'.
        :type action: str

        :return: Returns template path to log file.
        """
        self._validate_action(action)
        if action in ("preprocess_intervals", "annotate_gc"):
            name_pattern = "gcnv_{action}.{{library_kit}}".format(action=action)
            return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)
        elif action in (
            "filter_intervals",
            "contig_ploidy",
            "scatter_intervals",
            "merge_cohort_vcfs",
        ):
            name_pattern = "{{mapper}}.gcnv_{action}.{{library_kit}}".format(action=action)
            return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)
        elif action.startswith("call_cnvs"):
            name_pattern = "{{mapper}}.gcnv_{action}.{{library_kit}}.{{shard}}".format(
                action=action
            )
            return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)
        else:
            name_pattern = "{{mapper}}.gcnv_{action}.{{library_name}}".format(action=action)
            return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)

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
        )
        if action in high_resource_action_list:
            return self.resource_usage_dict.get("high_resource")
        else:
            return self.resource_usage_dict.get("default")

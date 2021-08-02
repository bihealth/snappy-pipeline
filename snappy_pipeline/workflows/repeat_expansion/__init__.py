# -*- coding: utf-8 -*
import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.repeat_expansion.annotate_expansionhunter import (
    AnnotateExpansionHunter,
)

#: Extensions of files to create as main payload - JSON.
EXT_JSON = (".json", ".json.md5")
#: Extensions of files to create as main payload - VCF.
EXT_VCF = (".vcf", ".vcf.md5")
#: Default configuration for the repeat_expansion step.
DEFAULT_CONFIG = r"""
# Default configuration repeat_expansion
step_config:
  repeat_expansion:
    # Repeat expansions definitions - used in ExpansionHunter call
    repeat_catalog: REQUIRED
    # Repeat expansions annotations, e.g., normality range - custom file
    repeat_annotation: REQUIRED
    # Path to the ngs_mapping step
    path_ngs_mapping: ../ngs_mapping
"""


class ExpansionHunterStepPart(BaseStepPart):
    """Repeat expansion analysis with Illumina::ExpansionHunter"""

    #: Step name.
    name = "expansionhunter"

    #: Valid actions.
    actions = ("run", "annotate")

    def get_input_files(self, action):
        """Return input function for ExpansionHunter rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns input function for ExpansionHunter rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate inputted action
        if action not in self.actions:
            valid_actions_str = ", ".join(self.actions)
            error_message = "Action '{action}' is not supported. Valid options: {options}".format(
                action=action, options=valid_actions_str
            )
            raise UnsupportedActionException(error_message)
        # Return requested function
        return getattr(self, "_get_input_files_{}".format(action))

    def get_output_files(self, action):
        """Return output function for ExpansionHunter rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns output function for ExpansionHunter rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate inputted action
        if action not in self.actions:
            valid_actions_str = ", ".join(self.actions)
            error_message = "Action '{action}' is not supported. Valid options: {options}".format(
                action=action, options=valid_actions_str
            )
            raise UnsupportedActionException(error_message)
        return getattr(self, "_get_output_files_{}".format(action))()

    def get_log_file(self, action):
        """Return log function for ExpansionHunter rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns log function for ExpansionHunter rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate inputted action
        if action not in self.actions:
            valid_actions_str = ", ".join(self.actions)
            error_message = "Action '{action}' is not supported. Valid options: {options}".format(
                action=action, options=valid_actions_str
            )
            raise UnsupportedActionException(error_message)
        return getattr(self, "_get_log_files_{}".format(action))()

    @listify
    def _get_input_files_run(self, wildcards):
        """Yield BAM files based on subworkflow `ngs_mapping` results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
        yield ngs_mapping(bam_tpl.format(**wildcards))

    @staticmethod
    @listify
    def _get_input_files_annotate(_wildcards):
        """Yield input files' pattern for rule `annotate` - based on ExpansionHunter call results.

        :param _wildcards: Snakemake rule wildcards (unused).
        :type _wildcards: snakemake.io.Wildcards
        """
        name_pattern = "{mapper}.expansionhunter.{library_name}"
        yield "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext="json"
        )

    @staticmethod
    @dictify
    def _get_output_files_run():
        """Yield output files' patterns for rule `run` - ExpansionHunter call."""
        # Initialise variables
        name_pattern = "{mapper}.expansionhunter.{library_name}"
        ext_dict = {"json": "json", "vcf": "vcf", "vcf_md5": "vcf.md5"}
        # Yield
        for key, ext in ext_dict.items():
            yield key, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    @staticmethod
    @dictify
    def _get_output_files_annotate():
        """Yield output files' patterns for rule `annotate`."""
        # Initialise variables
        name_pattern = "{mapper}.expansionhunter_annotated.{library_name}"
        ext_dict = {"json": "json", "json_md5": "json.md5"}
        # Yield
        for key, ext in ext_dict.items():
            yield key, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    @staticmethod
    def _get_log_files_run():
        """
        :return: Returns log file pattern for rule `run` - ExpansionHunter call.
        """
        name_pattern = "{mapper}.expansionhunter.{library_name}"
        return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)

    def annotate_results(self, _wildcards, sm_input, sm_output):
        """Annotate/Explain ExpansionHunter results.

        :param _wildcards: Snakemake wildcards associated with rule (unused).
        :type _wildcards: snakemake.io.Wildcards

        :param sm_input: Snakemake input associated with rule.
        :type sm_input: snakemake.io.Namedlist

        :param sm_output: Snakemake output associated with rule.
        :type sm_output: snakemake.io.Namedlist
        """
        # Absolute path from input and output
        input_path = os.path.join(os.getcwd(), str(sm_input))
        output_path = os.path.join(os.getcwd(), sm_output.json)
        # Annotate
        AnnotateExpansionHunter(
            eh_json=input_path,
            annotation_json=self.config["repeat_annotation"],
            output_path=output_path,
        ).run()


class RepeatExpansionWorkflow(BaseStep):
    """Perform germline repeat expansion analysis."""

    #: Workflow name
    name = "repeat_expansion"

    #: Sample sheet shortcut class
    sheet_shortcut_class = GermlineCaseSheet

    def __init__(
        self, workflow, config, cluster_config, config_lookup_paths, config_paths, workdir
    ):
        super().__init__(
            workflow,
            config,
            cluster_config,
            config_lookup_paths,
            config_paths,
            workdir,
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LinkOutStepPart, ExpansionHunterStepPart))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    @listify
    def _all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the germline repeat expansion analysis workflow."""
        # Initialise variable
        tools = ("expansionhunter",)
        # Yield the JSON annotated results files
        name_pattern = "{mapper}.{tool}_annotated.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            tool=tools,
            ext=EXT_JSON,
        )
        # Yield the VCF results files
        name_pattern = "{mapper}.{tool}.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            tool=tools,
            ext=EXT_VCF,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for donor in self._all_donors(include_background=False):
            yield from expand(tpl, donor=[donor], **kwargs)

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        # Requires path to ngs_mapping output, i.e., the BAM files
        self.ensure_w_config(
            config_keys=("step_config", "repeat_expansion", "path_ngs_mapping"),
            msg="Path to NGS mapping not configured but required for repeat expansion analysis.",
        )
        # Requires path to reference genome FASTA
        self.ensure_w_config(
            config_keys=("static_data_config", "reference", "path"),
            msg=(
                "Path to reference FASTA not configured but required "
                "for repeat expansion analysis."
            ),
        )

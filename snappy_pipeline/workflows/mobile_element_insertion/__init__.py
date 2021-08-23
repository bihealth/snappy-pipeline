# -*- coding: utf-8 -*
import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

#: Extensions of files to create as main payload.
EXT = ("_MEIs.txt", "_MEIs.txt.md5")


#: Default configuration for the mobile_element_insertion step.
DEFAULT_CONFIG = r"""
# Default configuration mobile_element_insertion
step_config:
  mobile_element_insertion:
    # Path to Scramble installation directory
    scramble_install_dir: REQUIRED
    # Path to the ngs_mapping step
    path_ngs_mapping: ../ngs_mapping
"""


class ScrambleStepPart(BaseStepPart):
    """Mobile element insertion detection using GeneDx::scramble"""

    #: Step name.
    name = "scramble"

    #: Valid actions.
    actions = ("cluster", "analysis")

    def get_input_files(self, action):
        """Return input function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns input function for scramble rule based on inputted action.

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
        """Return output function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns output function for scramble rule based on inputted action.

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

    @listify
    def _get_input_files_cluster(self, wildcards):
        """Yield BAM files based on subworkflow `ngs_mapping` results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bam_tpl = "output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam"
        yield ngs_mapping(bam_tpl.format(**wildcards))

    @staticmethod
    @listify
    def _get_input_files_analysis(wildcards):
        """Yield input files' pattern for rule `annotate` - based on scramble call results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        name_pattern = "{mapper}.scramble.{library_name}"
        base_name_out = "work/{name_pattern}/out/{name_pattern}_cluster.{ext}".format(
            name_pattern=name_pattern, ext="txt"
        )
        yield base_name_out.format(**wildcards)

    @staticmethod
    @dictify
    def _get_output_files_cluster():
        """Yield output files' patterns for scramble cluster call."""
        # Initialise variables
        name_pattern = "{mapper}.scramble.{library_name}"
        ext = "txt"
        # Yield
        yield ext, "work/{name_pattern}/out/{name_pattern}_cluster.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @staticmethod
    @dictify
    def _get_output_files_analysis():
        """Yield output files' patterns for scramble call."""
        # Initialise variables
        name_pattern = "{mapper}.scramble.{library_name}"
        ext_dict = {"txt": "_MEIs.txt", "txt_md5": "_MEIs.txt.md5"}
        # Yield
        for key, ext in ext_dict.items():
            yield key, "work/{name_pattern}/out/{name_pattern}{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    def get_log_file(self, action):
        """Return log function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns log function for scramble rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Initialise variable
        name_pattern = "{mapper}.scramble.{library_name}"

        # Validate inputted action
        if action not in self.actions:
            valid_actions_str = ", ".join(self.actions)
            error_message = "Action '{action}' is not supported. Valid options: {options}".format(
                action=action, options=valid_actions_str
            )
            raise UnsupportedActionException(error_message)
        if action == "annotate":
            name_pattern_annotated = "{mapper}.scramble_annotated.{library_name}"
            return "work/{name_pattern}/log/{name_pattern}.log".format(
                name_pattern=name_pattern_annotated
            )
        return "work/{name_pattern}/log/{name_pattern}_{action}.log".format(
            name_pattern=name_pattern, action=action
        )

    def get_params(self, action):
        """Get parameters.

        :param action: Action, i.e., step being performed.
        :type action: str

        :return: Returns method to get files required to run analysis part of scramble.
        """
        assert action == "analysis", "Parameters are only available for action 'analysis'."
        return self._get_analysis_parameters

    def _get_analysis_parameters(self, _wildcards):
        """Get parameters.

        :param _wildcards: Snakemake rule wildcards (unused).
        :type _wildcards: snakemake.io.Wildcards

        :return: Returns parameters required to run analysis part of scramble.
        """
        # Define required paths
        base_path = self.config["scramble_install_dir"]
        rscript_path = os.path.join(base_path, "SCRAMble.R")
        mei_refs_path = os.path.join(
            os.path.join(os.path.dirname(base_path), "resources"), "MEI_consensus_seqs.fa"
        )
        # Return dict with parameters
        return {"rscript": rscript_path, "mei_refs": mei_refs_path}


class MEIWorkflow(BaseStep):
    """Perform germline mobile element insertion detection."""

    #: Workflow name
    name = "mobile_element_insertion"

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
        self.register_sub_step_classes((LinkOutStepPart, ScrambleStepPart))
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
        """Get list of result files.

        :return: Return list of result files for the germline mobile element insertion
        detection workflow.
        """
        # Initialise variable
        tools = ("scramble",)
        name_pattern = "{mapper}.{tool}.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            tool=tools,
            ext=EXT,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for donor in self._all_donors(include_background=False):
            if donor.dna_ngs_library:  # ignores samples without DNA library
                yield from expand(tpl, donor=[donor], **kwargs)

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        # Requires path to ngs_mapping output, i.e., the BAM files
        self.ensure_w_config(
            config_keys=("step_config", "mobile_element_insertion", "path_ngs_mapping"),
            msg=(
                "Path to NGS mapping not configured but required for mobile "
                "element insertion detection."
            ),
        )
        # Requires path to ngs_mapping output, i.e., the BAM files
        self.ensure_w_config(
            config_keys=("step_config", "mobile_element_insertion", "scramble_install_dir"),
            msg=(
                "Path to scramble installation directory not configured but required for detection."
            ),
        )

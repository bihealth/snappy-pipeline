# -*- coding: utf-8 -*
"""Implementation of the ``targeted_seq_mei_calling`` step

The ``targeted_seq_mei_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs germline mobile element insertion identification (MEI).
The result are tabular files with mobile insertion characteristics (txt files).

==========
Stability
==========

This step is considered experimental, use it at your own discretion.

==========
Step Input
==========

MEI identification step uses Snakemake sub workflows for using the
result of the ``ngs_mapping`` step.

===========
Step Output
===========

For all samples, MEI identification will be performed on the primary DNA NGS libraries separately
for each configured read mapper and mobile element identification tool. The name of the primary DNA
NGS library will be used as an identification token in the output file.

For each read mapper, MEI tool, and sample, the following files will be generated:

- ``{mapper}.{mei_tool}.{lib_name}_MEIs.txt``
- ``{mapper}.{mei_tool}.{lib_name}_MEIs.txt.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.scramble.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.scramble.P001-N1-DNA1-WES1_MEIs.txt
    |       |-- bwa.scramble.P001-N1-DNA1-WES1_MEIs.txt.md5
    [...]


====================
Global Configuration
====================

Not applicable.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_targeted_seq_mei_calling.rst

==================================
Available MEI Identification Tools
==================================

The following germline MEI identification tool is currently available:

- ``"Scramble"``

=======
Reports
=======

Not applicable.

==================
Parallel Execution
==================

Not available.

"""
import os

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

#: Extensions of files to create as main payload.
EXT = ("_MEIs.txt", "_MEIs.txt.md5")


#: Default configuration for the targeted_seq_mei_calling step.
DEFAULT_CONFIG = r"""
# Default configuration
step_config:
  targeted_seq_mei_calling:
    # Full path to MEI reference file (fasta format)
    # if none is provided, it will use scramble's default
    mei_refs: null  # OPTIONAL
    # Minimum cluster size, depth of soft-clipped reads (set to Scramble default)
    n_cluster: 5  # OPTIONAL
    # Minimum MEI alignment score (set to Scramble default)
    mei_score: 50  # OPTIONAL
    # Minimum INDEL alignment score (set to Scramble default)
    indel_score: 80  # OPTIONAL
    # Minimum fraction of clipped length for calling polyA tail in MEIs (set to Scramble default)
    mei_polya_frac: 0.75  # OPTIONAL
    # Path to the ngs_mapping step
    path_ngs_mapping: ../ngs_mapping
"""


class ScrambleStepPart(BaseStepPart):
    """Mobile element insertion detection using GeneDx::scramble"""

    #: Step name
    name = "scramble"

    #: Class available actions
    actions = ("cluster", "analysis")

    def get_input_files(self, action):
        """Return input function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns input function for scramble rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action=action)
        return getattr(self, "_get_input_files_{}".format(action))

    def get_output_files(self, action):
        """Return output function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns output function for scramble rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action=action)
        return getattr(self, "_get_output_files_{}".format(action))()

    def get_log_file(self, action):
        """Return log function for scramble rules.

        :param action: Action (i.e., step) in the workflow.
        :type action: str

        :return: Returns log function for scramble rule based on inputted action.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action=action)
        # Set log
        name_pattern = "{mapper}.scramble.{library_name}"
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
        name_pattern = "{mapper}.scramble.{library_name}"
        ext = "txt"
        yield ext, "work/{name_pattern}/out/{name_pattern}_cluster.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @staticmethod
    @dictify
    def _get_output_files_analysis():
        """Yield output files' patterns for scramble call."""
        name_pattern = "{mapper}.scramble.{library_name}"
        ext_dict = {"txt": "_MEIs.txt", "txt_md5": "_MEIs.txt.md5"}
        for key, ext in ext_dict.items():
            yield key, "work/{name_pattern}/out/{name_pattern}{ext}".format(
                name_pattern=name_pattern, ext=ext
            )

    def _get_analysis_parameters(self, _wildcards):
        """Get parameters.

        :param _wildcards: Snakemake rule wildcards (unused).
        :type _wildcards: snakemake.io.Wildcards

        :return: Returns parameters required to run analysis part of scramble.
        """
        params = {
            "reference_genome": self.w_config["static_data_config"]["reference"]["path"],
            "mei_refs": self.config["mei_refs"],
            "n_cluster": self.config["n_cluster"],
            "mei_score": self.config["mei_score"],
            "indel_score": self.config["indel_score"],
            "mei_polya_frac": self.config["mei_polya_frac"],
        }
        return params

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="06:00:00",  # 6 hours
            memory=f"{8 * 1024}M",
        )


class MEIWorkflow(BaseStep):
    """Perform germline mobile element insertion detection."""

    #: Workflow name
    name = "targeted_seq_mei_calling"

    #: Sample sheet shortcut class
    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
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
            config_keys=("step_config", "targeted_seq_mei_calling", "path_ngs_mapping"),
            msg=(
                "Path to NGS mapping not configured but required for mobile "
                "element insertion detection."
            ),
        )
        # Requires reference FASTA to generate VCF files
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )

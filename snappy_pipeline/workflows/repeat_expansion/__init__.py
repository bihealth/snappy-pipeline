# -*- coding: utf-8 -*
"""Implementation of the ``repeat_analysis`` step

The ``repeat_analysis`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs repeat expansion analysis.  The result are variant files
(VCF) with the repeat expansions definitions, and associated annotations (JSON).

==========
Stability
==========

This step is considered experimental, use it at your own discretion.

==========
Step Input
==========

The repeat analysis step uses Snakemake sub workflows for using the result of the ``ngs_mapping``
step.

===========
Step Output
===========

For all samples, repeat analysis will be performed on the primary DNA NGS libraries separately for
each configured read mapper and repeat analysis tool. The name of the primary DNA NGS library will
be used as an identification token in the output file.

For each read mapper, repeat analysis tool, and sample, the following files will be generated:

- ``{mapper}.{repeat_tool}.{lib_name}.vcf``
- ``{mapper}.{repeat_tool}.{lib_name}.vcf.md5``
- ``{mapper}.{repeat_tool}_annotated.{lib_name}.json``
- ``{mapper}.{repeat_tool}_annotated.{lib_name}.json.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.expansionhunter.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.expansionhunter.P001-N1-DNA1-WES1.vcf
    |       |-- bwa.expansionhunter.P001-N1-DNA1-WES1.vcf.md5
    +-- bwa.expansionhunter_annotated.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.expansionhunter_annotated.P001-N1-DNA1-WES1.json
    |       |-- bwa.expansionhunter_annotated.P001-N1-DNA1-WES1.json.md5
    [...]

====================
Global Configuration
====================

Not applicable.

=====================
Default Configuration
=====================

The default configuration is as follows:

.. include:: DEFAULT_CONFIG_repeat_expansion.rst

===============================
Available Repeat Analysis Tools
===============================

The following germline repeat analysis tool is currently available:

- ``"ExpansionHunter"``


==================
Parallel Execution
==================

Not available.
"""

import os
from collections import OrderedDict

from biomedsheets.shortcuts import KEY_SEX, GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import UnsupportedActionException
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.repeat_expansion.annotate_expansionhunter import (
    AnnotateExpansionHunter,
)

from .model import RepeatExpansion as RepeatExpansionConfigModel

#: Extensions of files to create as main payload - JSON.
EXT_JSON = (".json", ".json.md5")
#: Extensions of files to create as main payload - VCF.
EXT_VCF = (".vcf", ".vcf.md5")
#: Default configuration for the repeat_expansion step.
DEFAULT_CONFIG = RepeatExpansionConfigModel.default_config_yaml_string()


class ExpansionHunterStepPart(BaseStepPart):
    """Repeat expansion analysis with Illumina::ExpansionHunter"""

    #: Step name
    name = "expansionhunter"

    #: Valid actions
    actions = ("run", "annotate")

    def __init__(self, *args, **kwargs):
        """Constructor."""
        super().__init__(*args, **kwargs)
        #: Build shortcut from library name to sex
        self.library_name_to_sex = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_sex.update(self._library_name_to_sex(sheet))

    @staticmethod
    def _library_name_to_sex(sheet):
        """Library name to sex.

        :param sheet: Sample sheet.
        :type sheet: biomedsheets.shortcuts.GermlineCaseSheet

        :return: Yields (library name, sex).
        """
        for donor in sheet.donors:
            sex = donor.extra_infos.get(KEY_SEX)
            for bio_sample in donor.bio_samples.values():
                for test_sample in bio_sample.test_samples.values():
                    for ngs_library in test_sample.ngs_libraries.values():
                        yield ngs_library.name, sex

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

    def _get_input_files_run(self, wildcards):
        """Yield BAM files based on subworkflow `ngs_mapping` results.

        :param wildcards: Snakemake rule wildcards.
        :type wildcards: snakemake.io.Wildcards
        """
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bam_tpl = ngs_mapping("output/{mapper}.{library_name}/out/{mapper}.{library_name}.bam")
        bam = bam_tpl.format(**wildcards)
        return {
            "bam": bam,
            "bai": bam + ".bai",
            "reference": self.w_config.static_data_config.reference.path,
            "repeat_catalog": self.config.repeat_catalog,
        }

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
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
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
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    @staticmethod
    def _get_log_files_run():
        """
        :return: Returns log file pattern for rule `run` - ExpansionHunter call.
        """
        name_pattern = "{mapper}.expansionhunter.{library_name}"
        return "work/{name_pattern}/log/{name_pattern}.log".format(name_pattern=name_pattern)

    def get_args(self, action):
        """Get parameters.

        :param action: Action, i.e., step being performed.
        :type action: str

        :return: Returns method to get donor's sex.
        """
        assert action == "run", "Parameters is only available for action 'run'."
        return self._get_donor_sex

    def _get_donor_sex(self, wildcards):
        """Get donor's sex.

        :param wildcards: Snakemake wildcards associated with rule (unused).
        :type wildcards: snakemake.io.Wildcards

        :return: Returns donor's sex as found in sample sheet: 'female', 'male' or 'unknown'.
        """
        return {"sex": self.library_name_to_sex[wildcards.library_name]}

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
            annotation_json=self.config.repeat_annotation,
            output_path=output_path,
        ).run()


class RepeatExpansionWorkflow(BaseStep):
    """Perform germline repeat expansion analysis."""

    #: Workflow name
    name = "repeat_expansion"

    #: Sample sheet shortcut class
    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=RepeatExpansionConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LinkOutStepPart, ExpansionHunterStepPart))
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

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
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            tool=tools,
            ext=EXT_JSON,
        )
        # Yield the VCF results files
        name_pattern = "{mapper}.{tool}.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            tool=tools,
            ext=EXT_VCF,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for donor in self._all_donors(include_background=False):
            if donor.dna_ngs_library:  # ignores samples without DNA library
                yield from expand(tpl, donor=[donor], **kwargs)

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        # Requires path to reference genome FASTA
        self.ensure_w_config(
            config_keys=("static_data_config", "reference", "path"),
            msg=(
                "Path to reference FASTA not configured but required for repeat expansion analysis."
            ),
        )

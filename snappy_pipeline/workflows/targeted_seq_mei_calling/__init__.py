# -*- coding: utf-8 -*
"""Implementation of the ``targeted_seq_mei_calling`` step

The ``targeted_seq_mei_calling`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and performs germline mobile element insertion (MEI) identification.
The result are VCF files with mobile insertions.

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

For each read mapper, MEI tool, and sample the following files will be generated:

- ``{mapper}.{mei_tool}.{lib_name}.vcf.gz``
- ``{mapper}.{mei_tool}.{lib_name}.vcf.gz.md5``

For example, it might look as follows for the example from above:

::

    output/
    +-- bwa.scramble.P001-N1-DNA1-WES1
    |   `-- out
    |       |-- bwa.scramble.P001-N1-DNA1-WES1.vcf.gz
    |       |-- bwa.scramble.P001-N1-DNA1-WES1.vcf.gz.md5
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

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import TargetedSeqMeiCalling as TargetedSeqMeiCallingConfigModel

#: Extensions of files to create as main payload.
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")


#: Default configuration for the targeted_seq_mei_calling step.
DEFAULT_CONFIG = TargetedSeqMeiCallingConfigModel.default_config_yaml_string()


class ScrambleStepPart(BaseStepPart):
    """Mobile element insertion detection using GeneDx::scramble"""

    #: Step name
    name = "scramble"

    #: Class available actions
    actions = ("cluster", "analysis")

    def __init__(self, parent):
        super().__init__(parent)

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

    def get_args(self, action):
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
        yield (
            ext,
            "work/{name_pattern}/out/{name_pattern}_cluster.{ext}".format(
                name_pattern=name_pattern, ext=ext
            ),
        )

    @staticmethod
    @dictify
    def _get_output_files_analysis():
        """Yield output files' patterns for scramble call."""
        name_pattern = "{mapper}.scramble.{library_name}"
        ext_dict = {
            "txt": "_MEIs.txt",
            "txt_md5": "_MEIs.txt.md5",
            "vcf": ".vcf",
            "vcf_gz": ".vcf.gz",
            "vcf_gz_md5": ".vcf.gz.md5",
            "vcf_tbi": ".vcf.gz.tbi",
            "vcf_tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in ext_dict.items():
            yield (
                key,
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                ),
            )

    def _get_analysis_parameters(self, _wildcards):
        """Get parameters.

        :param _wildcards: Snakemake rule wildcards (unused).
        :type _wildcards: snakemake.io.Wildcards

        :return: Returns parameters required to run analysis part of scramble.

        :raises InvalidConfiguration: if information provided in configuration isn't enough to run
        the analysis.
        """
        blast_ref_path = self.config.scramble.blast_ref
        try:
            if not os.path.isfile(blast_ref_path):
                raise InvalidConfiguration(
                    f"Provided path to reference genome ('blast_ref') "
                    f"is not a real file: {blast_ref_path}"
                )
        except TypeError as e:
            raise TypeError("Path to reference genome ('blast_ref') cannot be empty.") from e
        params = {
            "reference_genome": blast_ref_path,
            "mei_refs": self.config.scramble.mei_refs,
            "n_cluster": self.config.scramble.n_cluster,
            "mei_score": self.config.scramble.mei_score,
            "indel_score": self.config.scramble.indel_score,
            "mei_polya_frac": self.config.scramble.mei_polya_frac,
        }
        return params

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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


class MeiWorkflow(BaseStep):
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
            config_model_class=TargetedSeqMeiCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LinkOutStepPart, ScrambleStepPart))
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
        """Get list of result files.

        :return: Return list of result files for the germline mobile element insertion
        detection workflow.
        """
        # Initialise variable
        tools = ("scramble",)
        name_pattern = "{mapper}.{tool}.{donor.dna_ngs_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
            tool=tools,
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        for donor in self._all_donors(include_background=False):
            if donor.dna_ngs_library:  # ignores samples without DNA library
                yield from expand(tpl, donor=[donor], **kwargs)

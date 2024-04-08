# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_msi_calling`` step

This step allows for the detection of microsatellite instability for cancer samples from
whole genomes, exomes or large panels).  MANTIS starts from the aligned reads
(thus off ``ngs_mapping``) and generates a result file per tumor/normal pair.

==========
Step Input
==========

``somatic_msi_calling`` starts off the aligned reads, i.e. ``ngs_mapping``.

===========
Step Output
===========

Generally, the following links are generated to ``output/``.

.. note:: Tool-Specific Output

    As the only integrated tool is MANTIS at the moment, the output is very tailored to the result
    of this tool.  In the future, this section might contain "common" output and tool-specific
    output sub sections.

- ``mantis.{mapper}.{lib_name}-{lib_pk}/out/``
    - ``mantis.{mapper}.{lib_name}-{lib_pk}.results.txt``
    - ``mantis.{mapper}.{lib_name}-{lib_pk}.results.txt.status``

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_msi_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``mantis``

"""

import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: Default configuration for the somatic_msi_calling step
DEFAULT_CONFIG = r"""
# Default configuration somatic_msi_calling
step_config:
  somatic_msi_calling:
    tools: ['mantis']  # REQUIRED - available: 'mantis'
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    loci_bed: /fast/groups/cubi/projects/biotools/Mantis/appData/hg19/loci.bed  # REQUIRED
"""


class MantisStepPart(BaseStepPart):
    """Perform somatic microsatellite instability with MANTIS"""

    #: Step name
    name = "mantis"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/mantis.{{mapper}}.{{library_name}}/out/" "mantis.{{mapper}}.{{library_name}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.library_name]
        return pair.normal_sample.dna_ngs_library.name

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(
                    normal_library=self.get_normal_lib_name(wildcards), **wildcards
                )
            )
            tumor_base_path = (
                "output/{mapper}.{library_name}/out/" "{mapper}.{library_name}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        assert action == "run", "Unsupported actions"
        return input_function

    @dictify
    def get_output_files(self, action):
        assert action == "run", "Unsupported actions"
        exts = {"result": "results.txt", "status": "results.txt.status"}
        for key, ext in exts.items():
            yield (
                key,
                (
                    "work/mantis.{{mapper}}.{{library_name}}/out/"
                    "mantis.{{mapper}}.{{library_name}}_{sfx}"
                ).format(sfx=ext),
            )

    @staticmethod
    def get_log_file(action):
        """Return path to log file for the given action"""
        _ = action
        return "work/mantis.{mapper}.{library_name}/log/snakemake.mantis_run.log"

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=3,
            time="08:00:00",  # 2 hours
            memory=f"{30 * 1024 * 3}M",
        )


class SomaticMsiCallingWorkflow(BaseStep):
    """Perform somatic microsatellite instability analysis"""

    #: Step name
    name = "somatic_msi_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

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
        self.register_sub_step_classes((MantisStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the somatic targeted sequencing CNV calling step"""
        tool_actions = {"mantis": ("run",)}
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                for tool in self.config["tools"]:
                    for action in tool_actions[tool]:
                        try:
                            tpls = self.sub_steps[tool].get_output_files(action).values()
                        except AttributeError:
                            tpls = self.sub_steps[tool].get_output_files(action)
                        for tpl in tpls:
                            filenames = expand(
                                tpl,
                                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                                library_name=[sample_pair.tumor_sample.dna_ngs_library.name],
                            )
                            for f in filenames:
                                if ".tmp." not in f:
                                    yield f.replace("work/", "output/")

    def check_config(self):
        """Check that the necessary globalc onfiguration is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file not configured but required",
        )
        self.ensure_w_config(
            ("step_config", "somatic_msi_calling", "loci_bed"),
            "Path to bed file with microsatellite loci needed",
        )

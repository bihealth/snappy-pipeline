# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_msi_calling`` step

This step allows for the detection of microsatellite instability for cancer samples from
whole genomes, exomes or large panels).  MANTIS starts from the aligned reads
(thus off ``ngs_mapping``) and generates a result file per tumor/normal pair.

As MANTIS is not maintained anymore, the pipeline now supports only
`MANTIS2 <https://github.com/nh13/MANTIS2>`_.
The new version appears to be very silimar to the old one, both in terms of input & output files,
and in terms of requirements.

==========
Step Input
==========

``somatic_msi_calling`` starts off the aligned reads, i.e. ``ngs_mapping``.

===========
Step Output
===========

Generally, the following links are generated to ``output/``.

.. note:: Tool-Specific Output

    As the only integrated tool is MANTIS2 at the moment, the output is very tailored to the result
    of this tool.  In the future, this section might contain "common" output and tool-specific
    output sub sections.

- ``{mapper}.mantis_msi2.{lib_name}-{lib_pk}/out/``
    - ``{mapper}.mantis_msi2.{lib_name}-{lib_pk}.results.txt``
    - ``{mapper}.mantis_msi2.{lib_name}-{lib_pk}.results.txt.status``

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_msi_calling.rst

=====================================
Available Somatic Targeted CNV Caller
=====================================

- ``mantis_msi2``

"""

import os
import sys
from collections import OrderedDict
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

from .model import SomaticMsiCalling as SomaticMsiCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".results.txt", ".results.txt.status", ".results.txt.md5", ".results.txt.status.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("result", "status", "result_md5", "status_md5")

EXT_MATCHED = {
    "mantis_msi2": {
        "result": ".results.txt",
        "status": ".results.txt.status",
        "result_md5": ".results.txt.md5",
        "status_md5": ".results.txt.status.md5",
    },
}

#: Available somatic variant callers assuming matched samples.
MSI_CALLERS_MATCHED = ("mantis_msi2",)


#: Default configuration for the somatic_msi_calling step
DEFAULT_CONFIG = SomaticMsiCallingConfigModel.default_config_yaml_string()


class Mantis2StepPart(BaseStepPart):
    """Perform somatic microsatellite instability with MANTIS_msi2"""

    name = "mantis_msi2"

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{msi_caller}.{{tumor_library}}/out/"
            "{{mapper}}.{msi_caller}.{{tumor_library}}{ext}"
        )
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

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
                "output/{mapper}.{tumor_library}/out/{mapper}.{tumor_library}"
            ).format(**wildcards)
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
            }

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, msi_caller=[self.name], ext=EXT_VALUES))
        )

    def get_args(self, action: str) -> dict[str, Any]:
        # Validate action
        self._validate_action(action)
        return {
            "reference": self.parent.w_config.static_data_config.reference.path,
            "loci_bed": self.config.loci_bed,
        }

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)

        prefix = (
            "work/{{mapper}}.{msi_caller}.{{tumor_library}}/log/"
            "{{mapper}}.{msi_caller}.{{tumor_library}}"
        ).format(msi_caller=self.__class__.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=3,
            time="08:00:00",  # 8 hours
            memory=f"{30 * 1024 * 3}M",
        )


class SomaticMsiCallingWorkflow(BaseStep):
    """Perform somatic microsatellite instability analysis"""

    #: Step name
    name = "somatic_msi_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticMsiCallingConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((Mantis2StepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

    @listify
    def get_result_files(self):
        """Return list of result files for the MSI calling workflow"""
        name_pattern = "{mapper}.{msi_caller}.{tumor_library.name}"
        for msi_caller in set(self.config.tools) & set(MSI_CALLERS_MATCHED):
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                msi_caller=msi_caller,
                ext=EXT_MATCHED[msi_caller].values() if msi_caller in EXT_MATCHED else EXT_VALUES,
            )
            yield from self._yield_result_files_matched(
                os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
                mapper=self.w_config.step_config["ngs_mapping"].tools.dna,
                msi_caller=msi_caller,
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched msi callers such as
        mantis.
        """
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
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def check_config(self):
        """Check that the necessary globalc onfiguration is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file not configured but required",
        )

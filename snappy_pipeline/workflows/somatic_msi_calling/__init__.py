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

- ``{lib_name}-{lib_pk}/out/``
    - ``{lib_name}-{lib_pk}.results.txt``
    - ``{lib_name}-{lib_pk}.results.txt.status``

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

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments

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
        self.base_path_out = "work/{tumor_library}/out/{tumor_library}{ext}"
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
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_lib = self.get_normal_lib_name(wildcards)
            tumor_lib = wildcards.tumor_library
            normal: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=normal_lib
            )
            tumor: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=tumor_lib
            )
            return {
                "normal_bam": normal.bam,
                "normal_bai": normal.bai,
                "tumor_bam": tumor.bam,
                "tumor_bai": tumor.bai,
                "reference": self.w_config.static_data_config.reference.path,
                "loci_bed": self.config.loci_bed,
            }

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        return {
            key: self.base_path_out.replace("{ext}", ext) for key, ext in zip(EXT_NAMES, EXT_VALUES)
        }

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        # Validate action
        self._validate_action(action)

        prefix = "work/{tumor_library}/log/{tumor_library}"
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
            threads=1,  # Because of https://github.com/OSU-SRLab/MANTIS/issues/57
            runtime="24h",  # 24 hours
            mem=f"{30 * 1024 * 3}MB",
        )


class SomaticMsiCallingWorkflow(BaseStep):
    """Perform somatic microsatellite instability analysis"""

    #: Step name
    name = "somatic_msi_calling"
    consumes = {DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"msi"}))]

    config_model_class = SomaticMsiCallingConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local MSI calling output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"results": f"output/{lib}/out/{lib}.results.txt"}

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str | None = None,
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=(NgsMappingWorkflow,),
            task_name=task_name,
            **kwargs,
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((Mantis2StepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the MSI calling workflow"""
        msi_tool = str(self.config.tool)
        if msi_tool not in MSI_CALLERS_MATCHED:
            return
        name_pattern = "{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            ext=EXT_MATCHED[msi_tool].values() if msi_tool in EXT_MATCHED else EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
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
                    tpl,
                    tumor_library=[sample_pair.tumor_sample.dna_ngs_library],
                    **kwargs,
                )

    def check_config(self):
        """Check that the necessary globalc onfiguration is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA file not configured but required",
        )

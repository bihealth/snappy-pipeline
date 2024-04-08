# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_hla_loh_calling`` step

This step allows for the detection of loss of heterzygosity for cancer samples
from whole genomes, exomes or large panels).
LOHHLA starts from the aligned reads, germline HLA calls and optionally purity
and plodiy estimates of a sample.

==========
Step Input
==========

``somatic_hla_loh_calling`` starts off the aligned reads, i.e. ``ngs_mapping``,
HLA calls from ``hla_calling`` and results from
``somatic_purity_ploidy_estimate``

===========
Step Output
===========

A report file.
"""

import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: Default configuration for the somatic_msi_calling step
DEFAULT_CONFIG = r"""
# Default configuration somatic_hla_loh_calling
step_config:
  somatic_hla_loh_calling:
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    path_hla_typing: ../hla_typing  # REQUIRED
    path_somatic_purity_ploidy: ../somatic_purity_ploidy_estimate  # REQUIRED
"""


class LohhlaStepPart(BaseStepPart):
    """Perform LOHHLA analysis of tumor/normal WES/WGS pairs."""

    #: Step name
    name = "lohhla"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{{mapper}}.{{hla_caller}}.lohhla.{{tumor_library}}/out/"
            "{{mapper}}.{{hla_caller}}.lohhla.{{tumor_library}}{ext}"
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
            hla_typing = self.parent.sub_workflows["hla_typing"]
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
            hla = "output/optitype.{normal_library}/out/optitype.{normal_library}.txt".format(
                normal_library=self.get_normal_lib_name(wildcards)
            )

            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
                "tumor_bam": ngs_mapping(tumor_base_path + ".bam"),
                "tumor_bai": ngs_mapping(tumor_base_path + ".bam.bai"),
                "hla": hla_typing(hla),
            }

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        """Return output files from LOHHLA"""
        # Validate action
        self._validate_action(action)
        return {"done": expand(self.base_path_out, ext=".done")}

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = (
            "work/{mapper}.{hla_caller}.lohhla.{tumor_library}/log/"
            "{mapper}.{hla_caller}.lohhla.{tumor_library}"
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class SomaticHlaLohCallingWorkflow(BaseStep):
    """Perform somatic hla loh calling"""

    #: Workflow name
    name = "somatic_hla_loh_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
        self.register_sub_step_classes((LohhlaStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        self.register_sub_workflow("hla_typing", self.config["path_hla_typing"])

    @listify
    def get_result_files(self):
        """Return list of result files for the somatic hla loh calling workflow.

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        name_pattern = "{mapper}.optitype.lohhla.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
            ext=".done",
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
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

        This function returns the results from the matched somatic variant callers such as
        Mutect.
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
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_hla_loh_calling", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required.",
        )

# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_signatures`` step

The ``somatic_variant_signatures`` step takes as the input the results of the
``somatic_variant_calling`` step (bgzip-ed and indexed VCF files) and performs
deconstruction of signatures of the mutational processes at play.
The result it a data.frame tsv with the fraction of variants that each
signature explains as well as a plot.
"""


import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand
from collections import OrderedDict

from ..abstract import BaseStepPart, BaseStep, LinkOutStepPart
from ..ngs_mapping import NgsMappingWorkflow
from ..somatic_variant_calling import SomaticVariantCallingWorkflow, SOMATIC_VARIANT_CALLERS_MATCHED
from ...utils import listify, dictify

__author__ = "Clemens Messerschmidt"

# Default configuration variant_signatures
DEFAULT_CONFIG = r"""
step_config:
  somatic_variant_signatures:
    path_somatic_variant_calling: ../somatic_variant_calling   # REQUIRED
"""


class SignaturesStepPart(BaseStepPart):
    """Base class for signature classes"""

    def __init__(self, parent):
        super().__init__(parent)
        self.log_path = (
            "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/"
            "log/snakemake.tabulate_vcf.log"
        )
        # Build shortcut from cancer bio sample name to matched cancre sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor


class TabulateVariantsStepPart(SignaturesStepPart):
    """Tabulate mutation from VCF"""

    name = "tabulate_vcf"

    @dictify
    def get_input_files(self, action):
        """Return path to input file"""
        assert action == "run"
        tpl = (
            "output/{mapper}.{var_caller}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{tumor_library}"
        )
        KEY_EXT = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["somatic_variant_calling"]
        for key, ext in KEY_EXT.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files to tabulate vcf"""
        assert action == "run"
        yield "tsv", (
            "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/out/"
            "{mapper}.{var_caller}.tabulate_vcf.{tumor_library}.tsv"
        )

    def get_log_file(self, action):
        assert action == "run"
        return self.log_path

    def get_params(self, action):
        """Return arguments to pass down."""

        def params_function(wildcards):
            if wildcards.tumor_library not in self.donors:
                return {
                    "tumor_library": wildcards.tumor_library,
                    "normal_library": self.get_normal_lib_name(wildcards),
                }
            else:
                return {}

        return params_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["somatic_variant_signatures_tabulate_vcf"] = {
            "mem": 7 * 1024 * 2,
            "time": "01:00",
            "ntasks": 2,
        }


class DeconstructSigsStepPart(SignaturesStepPart):
    """Use deconstructSigs R package to identify signatures from tables"""

    name = "deconstruct_sigs"

    def __init__(self, parent):
        super().__init__(parent)
        self.log_path = (
            "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/"
            "log/snakemake.deconstruct_sigs.log"
        )

    @dictify
    def get_input_files(self, action):
        """Return input files to deconstruct signatures"""
        assert action == "run"
        yield "tsv", (
            "work/{mapper}.{var_caller}.tabulate_vcf.{tumor_library}/out/"
            "{mapper}.{var_caller}.tabulate_vcf.{tumor_library}.tsv"
        )

    @dictify
    def get_output_files(self, action):
        """Return output files to deconstruct signatures"""
        assert action == "run"
        yield "tsv", (
            "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/out/"
            "{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}.tsv"
        )
        yield "pdf", (
            "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/out/"
            "{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}.pdf"
        )

    def get_log_file(self, action):
        assert action == "run"
        return (
            "work/{mapper}.{var_caller}.deconstruct_sigs.{tumor_library}/"
            "log/snakemake.deconstruct_sigs.log"
        )

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["somatic_variant_signatures_deconstruct_sigs"] = {
            "mem": 7 * 1024 * 2,
            "time": "01:00",
            "ntasks": 2,
        }


class SomaticVariantSignaturesWorkflow(BaseStep):
    """Perform somatic variant signatures"""

    name = "somatic_variant_signatures"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

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
            (SomaticVariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (TabulateVariantsStepPart, DeconstructSigsStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow(
            "somatic_variant_calling", self.config["path_somatic_variant_calling"]
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.get("tools_ngs_mapping"):
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config.get("tools_somatic_variant_calling"):
            self.config["tools_somatic_variant_calling"] = self.w_config["step_config"][
                "somatic_variant_calling"
            ]["tools"]

    @listify
    def get_result_files(self):
        """Return list of result files for workflow"""
        callers = set(self.config["tools_somatic_variant_calling"])
        token = "{mapper}.{caller}.deconstruct_sigs.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", token, "out", token + ".tsv"),
            mapper=self.config["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
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
                        "normal or primary cancer library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library], **kwargs
                )

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_variant_signatures", "path_somatic_variant_calling"),
            ("Path to variant calling not configured but required for somatic variant signatures"),
        )

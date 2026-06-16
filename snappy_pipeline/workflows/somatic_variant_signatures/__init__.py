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
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import ResourceUsage
from snappy_pipeline.workflows.somatic_variant_calling.model import ExpectedSomaticVariants

from .model import SomaticVariantSignatures as SomaticVariantSignaturesConfigModel

__author__ = "Clemens Messerschmidt"


# Default configuration variant_signatures
DEFAULT_CONFIG = SomaticVariantSignaturesConfigModel.default_config_yaml_string()


class SignaturesStepPart(BaseStepPart):
    """Base class for signature classes"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)

        self.name_postfix = "{tumor_library}"

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

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        name_pattern = f"{self.name}." + self.name_postfix
        return os.path.join("work", name_pattern, "log", name_pattern + ".log")

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            runtime="1h",  # 1 hour
            mem=f"{7 * 1024 * 2}MB",
        )


class TabulateVariantsStepPart(SignaturesStepPart):
    """Tabulate mutation from VCF"""

    #: Step name
    name = "tabulate_vcf"

    @dictify
    def get_input_files(self, action):
        """Return path to input file"""
        # Validate action
        self._validate_action(action)
        name_pattern = self.name_postfix
        variants: ExpectedSomaticVariants = self.parent.get_upstream_paths(
            "somatic_variant", library_name=name_pattern
        )
        yield "vcf", variants.vcf
        yield "vcf_tbi", variants.vcf_tbi

    @dictify
    def get_output_files(self, action):
        """Return output files to tabulate vcf"""
        # Validate action
        self._validate_action(action)
        name_pattern = "tabulate_vcf." + self.name_postfix
        yield "tsv", os.path.join("work", name_pattern, "out", name_pattern + ".tsv")

    def get_args(self, action):
        """Return arguments to pass down."""
        # Validate action
        self._validate_action(action)

        def args_fn(wildcards):
            if wildcards.tumor_library not in self.donors:
                return {
                    "tumor_library": wildcards.tumor_library,
                    "normal_library": self.get_normal_lib_name(wildcards),
                }
            else:
                return {}

        return args_fn

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name


class DeconstructSigsStepPart(SignaturesStepPart):
    """Use deconstructSigs R package to identify signatures from tables"""

    #: Step name
    name = "deconstruct_sigs"

    def __init__(self, parent):
        super().__init__(parent)

    @dictify
    def get_input_files(self, action):
        """Return input files to deconstruct signatures"""
        # Validate action
        self._validate_action(action)
        name_pattern = "tabulate_vcf." + self.name_postfix
        yield "tsv", os.path.join("work", name_pattern, "out", name_pattern + ".tsv")

    @dictify
    def get_output_files(self, action):
        """Return output files to deconstruct signatures"""
        # Validate action
        self._validate_action(action)
        name_pattern = "deconstruct_sigs." + self.name_postfix
        yield "tsv", os.path.join("work", name_pattern, "out", name_pattern + ".tsv")
        yield "pdf", os.path.join("work", name_pattern, "out", name_pattern + ".pdf")


class SomaticVariantSignaturesWorkflow(BaseStep):
    """Perform somatic variant signatures"""

    #: Workflow name
    name = "somatic_variant_signatures"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"somatic", ("snv", "indel")})): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"signatures"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    config_model_class = SomaticVariantSignaturesConfigModel

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local signature output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"tsv": f"output/deconstruct_sigs.{lib}/out/deconstruct_sigs.{lib}.tsv"}

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
            previous_steps=(),
            task_name=task_name,
            **kwargs,
        )

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (TabulateVariantsStepPart, DeconstructSigsStepPart, LinkOutStepPart)
        )

    @listify
    def get_result_files(self):
        """Return list of result files for workflow"""
        name_pattern = "deconstruct_sigs.{tumor_library.name}"

        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + ".tsv")
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
                    tpl,
                    tumor_library=[sample_pair.tumor_sample.dna_ngs_library],
                    **kwargs,
                )

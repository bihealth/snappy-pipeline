# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_variant_signatures`` step

The ``somatic_variant_signatures`` step takes as the input the results of the
``somatic_variant_calling`` step (bgzip-ed and indexed VCF files) and performs
deconstruction of signatures of the mutational processes at play.
The result it a data.frame tsv with the fraction of variants that each
signature explains as well as a plot.
"""

from collections import OrderedDict
import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_annotation import (
    ANNOTATION_TOOLS,
    SomaticVariantAnnotationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_filtration import SomaticVariantFiltrationWorkflow

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
        self.config = parent.w_config.step_config["somatic_variant_signatures"]
        self.name_pattern = "{mapper}.{var_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                self.name_pattern += ".{anno_caller}.filtered"
            else:
                self.name_pattern += ".{anno_caller}.dkfz_bias_filter.eb_filter"
        self.name_pattern += "." + self.name + ".{tumor_library}"
        if self.config.is_filtered and len(self.config.filters) > 0:
            self.name_pattern += ".{filter}.{region}"
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
        return os.path.join("work", self.name_pattern, "log", "snakemake." + self.name + ".log")

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
            time="01:00:00",  # 1 hour
            memory=f"{7 * 1024 * 2}M",
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
        name_pattern = "{mapper}.{var_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                name_pattern += ".{anno_caller}.filtered"
            else:
                name_pattern += ".{anno_caller}.dkfz_bias_filter.eb_filter"
        name_pattern += ".{tumor_library}"
        if self.config.is_filtered and len(self.config.filters) > 0:
            name_pattern += ".{filter}.{region}"
        tpl = os.path.join("output", name_pattern, "out", name_pattern)
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["somatic_variant"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files to tabulate vcf"""
        # Validate action
        self._validate_action(action)
        yield "tsv", os.path.join("work", self.name_pattern, "out", self.name_pattern + ".tsv")

    def get_params(self, action):
        """Return arguments to pass down."""
        # Validate action
        self._validate_action(action)

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
        name_pattern = "{mapper}.{var_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                name_pattern += ".{anno_caller}.filtered"
            else:
                name_pattern += ".{anno_caller}.dkfz_bias_filter.eb_filter"
        name_pattern += ".tabulate_vcf.{tumor_library}"
        if self.config.is_filtered and len(self.config.filters) > 0:
            name_pattern += ".{filter}.{region}"
        yield "tsv", os.path.join("work", name_pattern, "out", name_pattern + ".tsv")

    @dictify
    def get_output_files(self, action):
        """Return output files to deconstruct signatures"""
        # Validate action
        self._validate_action(action)
        yield "tsv", os.path.join("work", self.name_pattern, "out", self.name_pattern + ".tsv")
        yield "pdf", os.path.join("work", self.name_pattern, "out", self.name_pattern + ".pdf")


class SomaticVariantSignaturesWorkflow(BaseStep):
    """Perform somatic variant signatures"""

    #: Workflow name
    name = "somatic_variant_signatures"

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
            config_model_class=SomaticVariantSignaturesConfigModel,
            previous_steps=(
                SomaticVariantCallingWorkflow,
                SomaticVariantAnnotationWorkflow,
                SomaticVariantFiltrationWorkflow,
                NgsMappingWorkflow,
            ),
        )
        # Register sub workflows
        config = self.w_config.step_config["somatic_variant_signatures"]
        sub_workflow = "somatic_variant_calling"
        if config.is_filtered:
            sub_workflow = "somatic_variant_filtration"
        self.register_sub_workflow(sub_workflow, config.path_somatic_variant, "somatic_variant")
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not config.tools_ngs_mapping:
            config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not config.tools_somatic_variant_calling:
            config.tools_somatic_variant_calling = self.w_config.step_config[
                "somatic_variant_calling"
            ].tools
        if not config.tools_somatic_variant_annotation:
            config.tools_somatic_variant_annotation = self.w_config.step_config[
                "somatic_variant_annotation"
            ].tools
        if config.is_filtered:
            if len(self.w_config.step_config["somatic_variant_filtration"].filter_list) > 0:
                config.filters = []
                config.filtered_regions = []
            else:
                if not config.filters:
                    config.filters = list(
                        self.w_config.step_config["somatic_variant_filtration"].filter_sets.keys()
                    )
                    config.filters.append("no_filter")
                if not config.filtered_regions:
                    config.filtered_regions = list(
                        self.w_config.step_config["somatic_variant_filtration"].exon_lists.keys()
                    )
                    config.filtered_regions.append("genome_wide")
        self.w_config.step_config["somatic_variant_signatures"] = config
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (TabulateVariantsStepPart, DeconstructSigsStepPart, LinkOutStepPart)
        )

    @listify
    def get_result_files(self):
        """Return list of result files for workflow"""
        config = self.w_config.step_config["somatic_variant_signatures"]
        name_pattern = "{mapper}.{caller}"
        if config.is_filtered:
            if len(config.filters) > 0:
                name_pattern += ".{anno_caller}.dkfz_bias_filter.eb_filter"
            else:
                name_pattern += ".{anno_caller}.filtered"
        name_pattern += ".deconstruct_sigs.{tumor_library.name}"
        if config.is_filtered and len(config.filters) > 0:
            name_pattern += ".{filter}.{region}"

        mappers = set(config.tools_ngs_mapping) & set(
            self.w_config.step_config["ngs_mapping"].tools.dna
        )
        assert len(mappers) > 0, "No valid mapper"
        callers = set(config.tools_somatic_variant_calling) & set(SOMATIC_VARIANT_CALLERS_MATCHED)
        assert len(callers) > 0, "No valid somatic variant caller"

        anno_callers = []
        filters = []
        regions = []
        if config.is_filtered:
            anno_callers = set(config.tools_somatic_variant_annotation) & set(ANNOTATION_TOOLS)
            assert len(anno_callers) > 0, "No valid somatic variant annotation tool"
            if len(config.filters) > 0:
                filters = list(
                    self.w_config.step_config["somatic_variant_filtration"].filter_sets.keys()
                )
                filters.append("no_filter")
                filters = set(filters) & set(config.filters)
                regions = list(
                    self.w_config.step_config["somatic_variant_filtration"].exon_lists.keys()
                )
                regions.append("genome_wide")
                regions = set(regions) & set(config.filtered_regions)

        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + ".tsv"),
            mapper=mappers,
            caller=callers,
            anno_caller=anno_callers,
            filter=filters,
            region=regions,
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
        if self.config.is_filtered:
            self.ensure_w_config(
                ("step_config", "somatic_variant_filtration"),
                "When is_filtered is set to True, "
                "the somatic_variant_filtration step must be configured",
            )

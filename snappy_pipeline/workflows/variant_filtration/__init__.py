# -*- coding: utf-8 -*-
"""Implementation of the ``variant_filtration`` step

This step takes annotated variants as the input from ``variant_annotation`` and performs various
filtration and postprocessing operations:

1. filter to high-confidence variants
    1. apply quality filter sets
    2. filter for consistency between different callers
2. filter to compatible mode of inheritance
3. filter by population/cohort frequency, remove polymorphisms
4. filter by region
5. filter by scores (e.g., conservation)
6. filter for het. comp. inheritance or keep all

# ::

#     1
#     stringent
#     loose

#     2
#     $qual.denovo
#     $qual.dom
#     $qual.rec_hom

#     3
#     $qual.denovo.denov_freq
#     $qual.dom.dom_freq
#     $qual.dom.rec_freq
#     $qual.rec_hom.rec_freq

#     4
#     $qual.denovo.denov_freq.$region
#     $qual.dom.dom_freq.$region
#     $qual.dom.rec_freq.$region
#     $qual.rec_hom.rec_freq.$region

#     5
#     $qual.denovo.denov_freq.$region.$scores
#     $qual.dom.dom_freq.$region.$scores
#     $qual.dom.rec_freq.$region.$scores
#     $qual.rec_hom.rec_freq.$region.$scores

#     6
#     $qual.denovo.denov_freq.$region.keep_all
#     $qual.dom.dom_freq.$region.keep_all
#     $qual.dom.rec_freq.$region.$scores.same_gene
#     $qual.dom.rec_freq.$region.$scores.same_tad
#     $qual.dom.rec_freq.$region.$scores.itv_500bp
#     $qual.rec_hom.rec_freq.$region.keep_all

================
Filtration Steps
================

The combinations of the filters is given in the configuration setting ``filter_combinations``
as dot-separated values, e.g., ``AA.BB.CC``.

==========
Step Input
==========

TODO

===========
Step Output
===========

TODO

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_filtration.rst

=======
Reports
=======

Currently, no reports are generated.
"""

# TODO: the implementation is super ugly and needs some refinement...

import os
import os.path
import sys
from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand, Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    InputFilesStepPartMixin,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow

from .model import VariantFiltration as VariantFiltrationConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = VariantFiltrationConfigModel.default_config_yaml_string()


class FiltersVariantsStepPartBase(BaseStepPart):
    """Base class for the different filters."""

    #: Step name
    name = None

    #: File name pattern
    name_pattern = None

    #: Class available actions
    actions = ("run",)

    #: Wildcard pattern name (must be set by derived class)
    filter_mode = None

    #: Model attribute name (must be set by derived class)
    filter_config = None

    def __init__(self, parent):
        super().__init__(parent)
        assert self.name_pattern is not None, "Set into class..."
        name_pattern = self.name_pattern
        self.base_path_out = os.path.join(
            "work", name_pattern, "out", name_pattern.replace(r",[^\.]+", "") + "{ext}"
        )
        self.path_log = os.path.join(
            "work", name_pattern, "out", name_pattern.replace(r",[^\.]+", "") + ".log"
        )

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
            time="1-00:00:00",  # 1 day
            memory=f"{int(3.75 * 1024 * 2)}M",
        )

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{ext}", ext)

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.path_log

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        def args_fn(wildcards: Wildcards) -> dict[str, Any]:
            assert self.filter_mode is not None, (
                f"'filter_mode' must be defined for sub-step '{self.name}"
            )
            params = {
                "index_library": wildcards.index_library,
                "filter_mode": getattr(wildcards, self.filter_mode),
            }
            if self.filter_config:
                params["filter_config"] = getattr(self.config, self.filter_config).model_dump(
                    by_alias=True
                )
            return params

        return args_fn


class FilterQualityStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_quality"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FiltersVariantsStepPartBase``
    prev_class = FiltersVariantsStepPartBase

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_quality
    filter_mode = "thresholds"

    #: Model name for filter_quality
    filter_config = "thresholds"

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        @dictify
        def input_function(wildcards):
            yield (
                "ped",
                os.path.realpath(
                    "work/write_pedigree.{index_library}/out/{index_library}.ped"
                ).format(**wildcards),
            )
            variant_annotation = self.parent.sub_workflows["variant_annotation"]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                output_path = (
                    "output/{mapper}.{caller}.jannovar_annotate_vcf.{index_library}/out/"
                    "{mapper}.{caller}.jannovar_annotate_vcf.{index_library}"
                ).format(**wildcards)
                yield key, variant_annotation(output_path) + ext

        return input_function


class FilterInheritanceStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_inheritance"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FilterQualityStepPart``
    prev_class = FilterQualityStepPart

    #: Include pedigree file flag (True)
    include_ped_file = True

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_inheritance
    filter_mode = "inheritance"


class FilterFrequencyStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_frequency"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FilterInheritanceStepPart``
    prev_class = FilterInheritanceStepPart

    #: Include pedigree file flag (True)
    include_ped_file = True

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_frequency
    filter_mode = "frequency"

    #: Model name for filter_frequency
    filter_config = "frequencies"


class FilterRegionsStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_regions"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FilterFrequencyStepPart``
    prev_class = FilterFrequencyStepPart

    #: Include pedigree file flag (True)
    include_ped_file = True

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_regions
    filter_mode = "regions"

    #: Model name for filter_regions
    filter_config = "region_beds"


class FilterScoresStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_scores"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FilterRegionsStepPart``
    prev_class = FilterRegionsStepPart

    #: Include pedigree file flag (True)
    include_ped_file = True

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_scores
    filter_mode = "scores"

    #: Model name for filter_scores
    filter_config = "score_thresholds"


class FilterHetCompStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    #: Step name
    name = "filter_het_comp"

    #: File name pattern
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}.{het_comp,[^\.]+}"
    )

    #: Pointer to the previous executed step, class ``FilterScoresStepPart``
    prev_class = FilterScoresStepPart

    #: Include pedigree file flag (True)
    include_ped_file = True

    #: Types of output files by extension
    ext_names = EXT_NAMES

    #: Output file extensions
    ext_values = EXT_VALUES

    #: Wildcards name for filter_het_comp
    filter_mode = "het_comp"

    #: Model name for filter_scores
    filter_config = "region_beds"


class VariantFiltrationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    #: Workflow name
    name = "variant_filtration"

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

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
            config_model_class=VariantFiltrationConfigModel,
            previous_steps=(VariantAnnotationWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                FilterQualityStepPart,
                FilterInheritanceStepPart,
                FilterFrequencyStepPart,
                FilterRegionsStepPart,
                FilterScoresStepPart,
                FilterHetCompStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("variant_annotation", self.config.path_variant_annotation)
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not self.config.tools_variant_calling:
            self.config.tools_variant_calling = self.w_config.step_config["variant_calling"].tools

    @listify
    def get_result_files(self):
        """Return list of result files for the variant filtration workflow."""
        # Generate output paths without extracting individuals.
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library.name}.{filters}"
        )
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config.tools_ngs_mapping,
            caller=self.config.tools_variant_calling,
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(msg.format(pedigree), file=sys.stderr)
                    continue
                yield from expand(
                    tpl,
                    index_library=[pedigree.index.dna_ngs_library],
                    filters=self.config.filter_combinations,
                    **kwargs,
                )

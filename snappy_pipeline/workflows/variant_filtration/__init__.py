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

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    InputFilesStepPartMixin,
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
step_config:
  variant_filtration:
    path_variant_annotation: ../variant_annotation
    tools_ngs_mapping: null      # defaults to ngs_mapping tool
    tools_variant_calling: null  # defaults to variant_annotation tool
    thresholds:                  # quality filter sets, "keep_all" implicitely defined
      conservative:
        min_gq: 40
        min_dp_het: 10
        min_dp_hom: 5
        include_expressions:
        - 'MEDGEN_COHORT_INCONSISTENT_AC=0'
      relaxed:
        min_gq: 20
        min_dp_het: 6
        min_dp_hom: 3
        include_expressions:
        - 'MEDGEN_COHORT_INCONSISTENT_AC=0'
    frequencies:             # values to use for frequency filtration
      af_dominant: 0.001     # AF (allele frequency) values
      af_recessive: 0.01
      ac_dominant: 3         # AC (allele count in gnomAD) values
    region_beds:             # regions to filter to, "whole_genome" implicitely defined
      all_tads: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed
      all_genes: /fast/projects/medgen_genomes/static_data/GRCh37/gene_bed/ENSEMBL_v75.bed.gz
      limb_tads: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      lifted_enhancers: /fast/projects/medgen_genomes/static_data/GRCh37/all_but_onlyMB.bed
      vista_enhancers: /fast/projects/medgen_genomes/static_data/GRCh37/vista_limb_enhancers.bed
    score_thresholds:        # thresholds on scores to filter to, "all_scores" implictely defined
      coding:
        require_coding: true
        require_gerpp_gt2: false
        min_cadd: null
      conservative:  # unused; TODO: rename?
        require_coding: false
        require_gerpp_gt2: false
        min_cadd: 0
      conserved:  # TODO: rename?
        require_coding: false
        require_gerpp_gt2: true
        min_cadd: null
    filter_combinations: # dot-separated {thresholds}.{inherit}.{freq}.{region}.{score}.{het_comp}
    - conservative.de_novo.dominant_freq.lifted_enhancers.all_scores.passthrough
    - conservative.de_novo.dominant_freq.lifted_enhancers.conserved.passthrough
    - conservative.de_novo.dominant_freq.limb_tads.all_scores.passthrough
    - conservative.de_novo.dominant_freq.limb_tads.coding.passthrough
    - conservative.de_novo.dominant_freq.limb_tads.conserved.passthrough
    - conservative.de_novo.dominant_freq.vista_enhancers.all_scores.passthrough
    - conservative.de_novo.dominant_freq.vista_enhancers.conserved.passthrough
    - conservative.de_novo.dominant_freq.whole_genome.all_scores.passthrough
    - conservative.de_novo.dominant_freq.whole_genome.coding.passthrough
    - conservative.de_novo.dominant_freq.whole_genome.conserved.passthrough
    - conservative.dominant.dominant_freq.lifted_enhancers.all_scores.passthrough
    - conservative.dominant.dominant_freq.lifted_enhancers.conserved.passthrough
    - conservative.dominant.dominant_freq.limb_tads.all_scores.passthrough
    - conservative.dominant.dominant_freq.limb_tads.coding.passthrough
    - conservative.dominant.dominant_freq.limb_tads.conserved.passthrough
    - conservative.dominant.dominant_freq.vista_enhancers.all_scores.passthrough
    - conservative.dominant.dominant_freq.vista_enhancers.conserved.passthrough
    - conservative.dominant.dominant_freq.whole_genome.all_scores.passthrough
    - conservative.dominant.dominant_freq.whole_genome.coding.passthrough
    - conservative.dominant.dominant_freq.whole_genome.conserved.passthrough
    - conservative.dominant.recessive_freq.lifted_enhancers.all_scores.intervals500
    - conservative.dominant.recessive_freq.lifted_enhancers.conserved.intervals500
    - conservative.dominant.recessive_freq.lifted_enhancers.conserved.tads
    - conservative.dominant.recessive_freq.limb_tads.all_scores.intervals500
    - conservative.dominant.recessive_freq.limb_tads.coding.gene
    - conservative.dominant.recessive_freq.limb_tads.conserved.intervals500
    - conservative.dominant.recessive_freq.limb_tads.conserved.tads
    - conservative.dominant.recessive_freq.vista_enhancers.all_scores.intervals500
    - conservative.dominant.recessive_freq.vista_enhancers.conserved.intervals500
    - conservative.dominant.recessive_freq.vista_enhancers.conserved.tads
    - conservative.dominant.recessive_freq.whole_genome.all_scores.intervals500
    - conservative.dominant.recessive_freq.whole_genome.coding.gene
    - conservative.dominant.recessive_freq.whole_genome.conserved.intervals500
    - conservative.dominant.recessive_freq.whole_genome.conserved.tads
    - conservative.recessive_hom.recessive_freq.lifted_enhancers.all_scores.passthrough
    - conservative.recessive_hom.recessive_freq.lifted_enhancers.conserved.passthrough
    - conservative.recessive_hom.recessive_freq.limb_tads.all_scores.passthrough
    - conservative.recessive_hom.recessive_freq.limb_tads.coding.passthrough
    - conservative.recessive_hom.recessive_freq.limb_tads.conserved.passthrough
    - conservative.recessive_hom.recessive_freq.vista_enhancers.all_scores.passthrough
    - conservative.recessive_hom.recessive_freq.vista_enhancers.conserved.passthrough
    - conservative.recessive_hom.recessive_freq.whole_genome.all_scores.passthrough
    - conservative.recessive_hom.recessive_freq.whole_genome.coding.passthrough
    - conservative.recessive_hom.recessive_freq.whole_genome.conserved.passthrough
    # The following are for input to variant_combination.
    - conservative.dominant.recessive_freq.whole_genome.coding.passthrough
    - conservative.dominant.recessive_freq.whole_genome.conserved.passthrough
"""


class FiltersVariantsStepPartBase(BaseStepPart):
    """Base class for the different filters."""

    #: Name of the step (e.g., for rule names)
    name = None
    #: Token to use in file name
    name_pattern = None

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

    def update_cluster_config(self, cluster_config):
        cluster_config["variant_filtration_{}_run".format(self.name)] = {
            "mem": int(3.75 * 1024 * 2),
            "time": "01:00:00",
            "ntasks": 2,
        }

    @dictify
    def get_output_files(self, action):
        assert action == "run"
        for key, ext in zip(EXT_NAMES, EXT_VALUES):
            yield key, self.base_path_out.replace("{ext}", ext)

    def get_log_file(self, action):
        assert action == "run"
        return self.path_log


class FilterQualityStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_quality"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}"
    )
    prev_class = FiltersVariantsStepPartBase
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            yield "ped", os.path.realpath(
                "work/write_pedigree.{index_library}/out/{index_library}.ped"
            ).format(**wildcards)
            variant_annotation = self.parent.sub_workflows["variant_annotation"]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                output_path = (
                    "output/{mapper}.{caller}.jannovar_annotate_vcf.{index_library}/out/"
                    "{mapper}.{caller}.jannovar_annotate_vcf.{index_library}"
                ).format(**wildcards)
                yield key, variant_annotation(output_path) + ext

        assert action == "run", "Unsupported actions"
        return input_function


class FilterInheritanceStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_inheritance"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}"
    )
    prev_class = FilterQualityStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class FilterFrequencyStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_frequency"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}"
    )
    prev_class = FilterInheritanceStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class FilterRegionsStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_regions"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}"
    )
    prev_class = FilterFrequencyStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class FilterScoresStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_scores"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}"
    )
    prev_class = FilterRegionsStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class FilterHetCompStepPart(InputFilesStepPartMixin, FiltersVariantsStepPartBase):
    """Apply the configured filters."""

    name = "filter_het_comp"
    name_pattern = (
        r"{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{frequency,[^\.]+}.{regions,[^\.]+}."
        r"{scores,[^\.]+}.{het_comp,[^\.]+}"
    )
    prev_class = FilterScoresStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class VariantFiltrationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "variant_filtration"
    sheet_shortcut_class = GermlineCaseSheet

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
            (VariantAnnotationWorkflow, NgsMappingWorkflow),
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
        self.register_sub_workflow("variant_annotation", self.config["path_variant_annotation"])
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]

    @listify
    def get_result_files(self):
        """Return list of result files for the variant filtration workflow."""
        # Generate output paths without extracting individuals.
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_vcf.filtered.{index_library.name}.{filters}"
        )
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_variant_calling"],
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
                    filters=self.config["filter_combinations"],
                    **kwargs
                )

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_filtration", "path_variant_annotation"),
            ("Path to variant_annotation not configured but required for variant_filtration"),
        )

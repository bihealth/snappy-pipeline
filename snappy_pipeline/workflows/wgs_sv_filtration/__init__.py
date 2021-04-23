# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_filtration`` step

This step takes annotated SVs as the input from ``wgs_sv_annotation`` and performs post-filtration
operations.

1. apply quality filter sets,
2. apply filters on the mode of inheritance, and
3. filter down the variants to specific regions.

Note that in contrast to the ``wgs_sv_annotation`` step, variants are hard-filtered in this step.
This means that variants are actually removed from the variant call file instead of just
flagging them not passing some filter criteria.

================
Filtration Steps
================

The combinations of the filters is given in the configuration setting ``filter_combinations``
as dot-separated values, e.g., ``no_filter.all.whole_genome``.

---------------
Quality Filters
---------------

The actual thresholds are configured with the ``thresholds`` settings.

The following filter sets are available

- ``no_filter``, no variant is removed based on quality
- ``conservative``, conservative requirements on evidence supporting a variant, suitable
  for the non-*de novo* case.
- ``conservative_de_novo`` conservative requiremnts on evidences supporting a variant, suitable
  for the *de novo* case.

-------------------
Inheritance Filters
-------------------

The following inheritance filter sets are defined.

.. note:: Currently, the inheritance filters only work properly for trios!

    The filters give sensible results for singletons.  However, in case of missing parents, the
    filters fall back to sensible defaults.  Additional information such as siblings or more
    relatives are currently ignored.

- ``all``, no variant is removed based on inheritance.
- ``dominant``, dominant mode of inheritance.
- ``de_novo``, dominant mode of inheritance.
- ``recessive_hc``, recessive mode of inheritance, heterozygous composite case.
- ``recessive_hom``, recessive mode of inheritance, homozygous case.

--------------
Region Filters
--------------

In the ``regions_bed``, a key/mapping can be configured for specifying name/BED path pairs of
defined region sets.

Further, the ``whole_genome`` region is defined and leads to no filtration of variants based
on genomic region.

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

.. include:: DEFAULT_CONFIG_wgs_sv_filtration.rst

=======
Reports
=======

Currently, no reports are generated.
"""

import os
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
from snappy_pipeline.workflows.wgs_sv_annotation import WgsSvAnnotationWorkflow
from snappy_pipeline.workflows.wgs_sv_calling import WgsSvCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration for wgs_sv_filtration
step_config:
  wgs_sv_filtration:
    path_wgs_sv_annotation: ../wgs_sv_annotation
    # Regions.
    #
    # The "all_tads" entry has to be present for the het. comp. combination of variants.
    region_beds:                    # regions to filter to, "whole_genome" is implicitely defined
      all_tads: /fast/projects/medgen_genomes/static_data/GRCh37/hESC_hg19_allTads.bed  # REQUIRED
      limb_tads: /fast/projects/medgen_genomes/static_data/GRCh37/newlimb_tads.bed
      lifted_enhancers: /fast/projects/medgen_genomes/static_data/GRCh37/all_but_onlyMB.bed
      vista_enhancers: /fast/projects/medgen_genomes/static_data/GRCh37/vista_limb_enhancers.bed
    tools_ngs_mapping: False        # defaults to those in NGS mapping
    tools_wgs_sv_calling: False     # defaults to those in WGS SV calling
    # Supported modes of inheritance:
    # - de_novo
    # - dominant
    # - dominant_rp (reduced penetrance)
    # - recessive_hc (heterozygous composite)
    # - recessive_hom (homozygous)
    thresholds:                     # quality filter set, 'no_filter' is implicit
      conservative:
        min_pe_aaf_index: 0.3
        min_sr_aaf_index: 0.3
        min_pe_aaf_parent: 0
        min_sr_aaf_parent: 0.3
      conservative_de_novo:
        min_pe_aaf_index: 0.3
        min_sr_aaf_index: 0.3
        min_pe_aaf_parent: 0
        min_sr_aaf_parent: 0
    # TODO: Currently, all combinations are used, not necessary.
    filter_combinations:  # dot-separated {thresholds}.{inheritance}.{regions}
    - conservative_de_novo.de_novo.lifted_enhancers
    - conservative_de_novo.de_novo.limb_tads
    - conservative_de_novo.de_novo.vista_enhancers
    - conservative_de_novo.de_novo.whole_genome
    - conservative.dominant.lifted_enhancers
    - conservative.dominant.limb_tads
    - conservative.dominant_rp.lifted_enhancers
    - conservative.dominant_rp.limb_tads
    - conservative.dominant_rp.vista_enhancers
    - conservative.dominant_rp.whole_genome
    - conservative.dominant.vista_enhancers
    - conservative.dominant.whole_genome
    - conservative.recessive_hc.lifted_enhancers
    - conservative.recessive_hc.limb_tads
    - conservative.recessive_hc.vista_enhancers
    - conservative.recessive_hc.whole_genome
    - conservative.recessive_hom.lifted_enhancers
    - conservative.recessive_hom.limb_tads
    - conservative.recessive_hom.vista_enhancers
    - conservative.recessive_hom.whole_genome
"""


class FiltersWgsSvStepPartBase(BaseStepPart):
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
        cluster_config["wgs_sv_filtration_{}_run".format(self.name)] = {
            "mem": int(3.75 * 1024 * 2),
            "time": "01:00",
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


class FilterQualityStepPart(InputFilesStepPartMixin, FiltersWgsSvStepPartBase):
    """Apply the configured filters."""

    name = "filter_quality"
    name_pattern = (
        r"{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}.{thresholds,[^\.]+}"
    )
    prev_class = FiltersWgsSvStepPartBase
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES

    def get_input_files(self, action):
        @dictify
        def input_function(wildcards):
            yield "ped", "work/write_pedigree.{index_library}/out/{index_library}.ped".format(
                **wildcards
            )
            wgs_sv_annotation = self.parent.sub_workflows["wgs_sv_annotation"]
            for key, ext in zip(EXT_NAMES, EXT_VALUES):
                output_path = (
                    "output/{mapper}.{caller}.annotated.{index_library}/out/"
                    "{mapper}.{caller}.annotated.{index_library}"
                ).format(**wildcards)
                yield key, wgs_sv_annotation(output_path) + ext

        assert action == "run", "Unsupported actions"
        return input_function


class FilterInheritanceStepPart(InputFilesStepPartMixin, FiltersWgsSvStepPartBase):
    """Apply the configured filters."""

    name = "filter_inheritance"
    name_pattern = (
        r"{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}"
    )
    prev_class = FilterQualityStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class FilterRegionsStepPart(InputFilesStepPartMixin, FiltersWgsSvStepPartBase):
    """Apply the configured filters."""

    name = "filter_regions"
    name_pattern = (
        r"{mapper}.{caller}.annotated.filtered.{index_library,[^\.]+}."
        r"{thresholds,[^\.]+}.{inheritance,[^\.]+}.{regions,[^\.]+}"
    )
    prev_class = FilterInheritanceStepPart
    include_ped_file = True
    ext_names = EXT_NAMES
    ext_values = EXT_VALUES


class WgsSvFiltrationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "wgs_sv_filtration"
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
            (WgsSvAnnotationWorkflow, WgsSvCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                FilterQualityStepPart,
                FilterInheritanceStepPart,
                FilterRegionsStepPart,
                LinkOutStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("wgs_sv_annotation", self.config["path_wgs_sv_annotation"])
        # Copy over "tools" setting from wgs_sv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_wgs_sv_calling"]:
            self.config["tools_wgs_sv_calling"] = self.w_config["step_config"]["wgs_sv_calling"][
                "tools"
            ]
        # Set implicitely defined config.
        self.set_default_config()

    @listify
    def get_result_files(self):
        """Return list of result files for the variant filtration workflow."""
        # Generate output paths without extracting individuals.
        name_pattern = "{mapper}.{caller}.annotated.filtered.{index_library.name}.{filters}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_sv_calling"]["dna"],
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
        """Check that the path to the WGS SV annotation step is present"""
        self.ensure_w_config(
            ("step_config", "wgs_sv_filtration", "path_wgs_sv_annotation"),
            ("Path to wgs_sv_annotation not configured but must be for " "wgs_sv_filtration"),
        )

    def set_default_config(self):
        """Set default configuration."""
        self.w_config["step_config"]["wgs_sv_filtration"]["thresholds"]["no_filter"] = {
            "min_pe_aaf_index": 0,
            "min_sr_aaf_index": 0,
            "min_pe_aaf_parent": 0,
            "min_sr_aaf_parent": 0,
        }

# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_annotation`` step

The ``wgs_sv_annotation`` step takes as the input the results of the ``wgs_sv_calling`` step
(called germline SVs) and ``variant_calling`` (called small germline variants) and performs
annotation and filtration of the structural variants.

Such filters include:

- quality filter for removing calls with low support,
- inheritance compatibility filter that checks for compatibility of inheritance (only works
  for trios currently),
- various filters related to lower the false discovery rate in rare/de novo variant calling,
  e.g., counting number of affected individuals in cohort outside the index' family.

.. note::

    Status: not implemented yet

==========
Step Input
==========

The variant annotation step uses the output of the following CUBI pipeline steps:

- ``wgs_sv_calling``
- ``variant_annotation``

===========
Step Output
===========

.. note: TODO

====================
Global Configuration
====================

.. note: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_wgs_sv_annotation.rst

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
    LinkOutStepPart,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Default configuration of the wgs_sv_filtration step
DEFAULT_CONFIG = r"""
# Default configuration wgs_sv_annotation
step_config:
  wgs_sv_annotation:
    path_ngs_mapping: ../ngs_mapping
    path_variant_calling: ../variant_calling
    path_wgs_sv_calling: ../wgs_sv_calling
    tool_ngs_mapping_variant_calling: bwa
    tool_variant_calling: gatk_hc
    tools_ngs_mapping:
    - bwa
    tools_wgs_sv_calling:
    - delly2
    path_alu_bed: ''
    path_db_bed: ''
"""


class VcfSvFilterStepPart(BaseStepPart):
    """Annotate VCF using wgs_sv_filter.py script."""

    name = "vcf_sv_filter"

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/.done"
        self.log_path = (
            "work/{mapper}.{caller}.annotated.{index_ngs_library}/"
            "log/snakemake.wgs_sv_filter.log"
        )

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        assert action == "run"

        @dictify
        def input_function(wildcards):
            if wildcards.caller == "popdel":
                name_pattern = "internal.concat_calls"
                key_ext = {"bcf": ".vcf.gz", "csi": ".vcf.gz.tbi"}
            else:  # delly2
                name_pattern = "merge_genotypes"
                key_ext = {"bcf": ".bcf", "csi": ".bcf.csi"}

            tpl = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
            yield "ped", tpl.format(**wildcards)
            # TODO: make work for non-Delly?
            tpl = (
                "work/{mapper}.{caller}.{index_ngs_library}/out/"
                "{mapper}.{caller}.{index_ngs_library}"
            )
            # SVs
            wgs_sv_calling = self.parent.sub_workflows["wgs_sv_calling"]
            for key, ext in key_ext.items():
                yield "sv_" + key, wgs_sv_calling(
                    tpl.replace("{index_ngs_library}", name_pattern) + ext
                ).format(**wildcards)

            # Small variants
            key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
            variant_calling = self.parent.sub_workflows["variant_calling"]
            for key, ext in key_ext.items():
                path = tpl.replace("{caller}", self.parent.config["tool_variant_calling"])
                yield "var_" + key, variant_calling(path + ext).format(**wildcards)

        return input_function

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        assert action == "run"
        prefix = (
            "work/{mapper}.{caller}.annotated.{index_ngs_library}/out/"
            "{mapper}.{caller}.annotated.{index_ngs_library}"
        )
        key_ext = {
            "vcf": ".vcf.gz",
            "tbi": ".vcf.gz.tbi",
            "vcf_md5": ".vcf.gz.md5",
            "tbi_md5": ".vcf.gz.tbi.md5",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    def get_log_file(self, action):
        assert action == "run"
        return self.log_path

    @classmethod
    def update_cluster_config(cls, cluster_config):
        """Update cluster configuration with resource requirements"""
        cluster_config["wgs_sv_annotation_wgs_sv_filter"] = {
            "mem": 5 * 1024 * 2,
            "time": "100:00",
            "ntasks": 2,
        }


class WgsSvAnnotationWorkflow(BaseStep):
    """Perform germline WGS SV annotation"""

    name = "wgs_sv_annotation"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
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
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, VcfSvFilterStepPart, LinkOutStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])
        self.register_sub_workflow("wgs_sv_calling", self.config["path_wgs_sv_calling"])
        # Copy over "tools" setting from wgs_sv_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_wgs_sv_calling"]:
            self.config["tools_wgs_sv_calling"] = self.w_config["step_config"]["wgs_sv_calling"][
                "tools"
            ]

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "{mapper}.{caller}.annotated.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            caller=self.config["tools_wgs_sv_calling"],
            ext=EXT_VALUES,
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(
                        msg.format(  # pragma: no cover
                            list(sorted(d.name for d in pedigree.donors))
                        ),
                        file=sys.stderr,
                    )
                    continue  # pragma: no cover
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "wgs_sv_annotation", "path_variant_calling"),
            (" not configured but required for WGS SV annotation"),
        )
        self.ensure_w_config(
            ("step_config", "wgs_sv_annotation", "tool_ngs_mapping_variant_calling"),
            (
                "NGS mapping tool (of variant calling) not configured but required for "
                "WGS SV annotation"
            ),
        )
        self.ensure_w_config(
            ("step_config", "wgs_sv_annotation", "tool_variant_calling"),
            (
                "NGS variant caller (of variant calling) not configured but required "
                "for WGS SV annotation"
            ),
        )
        self.ensure_w_config(
            ("step_config", "wgs_sv_annotation", "tools_ngs_mapping"),
            ("NGS mapping tools not configured but required for WGS SV annotation"),
        )
        self.ensure_w_config(
            ("step_config", "wgs_sv_annotation", "tools_wgs_sv_calling"),
            ("WGS SV calling tools not configured but required for WGS SV annotation"),
        )

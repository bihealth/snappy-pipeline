# -*- coding: utf-8 -*-
"""Implementation of the germline ``variant_checking`` step

The ``variant_checking`` step takes as the input the results of the ``variant_annotation``
step.  It then executes various tools computing statistics on the result files and consistency
checks with the pedigrees.

==========
Step Input
==========

The variant calling step uses Snakemake sub workflows for using the result of the
``variant_annotation`` step.

===========
Step Output
===========

.. note:: TODO

====================
Global Configuration
====================

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_checking.rst

==========================
Available Variant Checkers
==========================

The following variant checkers integrated:

- ``"bcftools_stats"`` -- call ``bcftools stats`` for various statistics
- ``"peddy"`` -- check variants against a PED file

=======
Reports
=======

Currently, no reports are generated.
"""

import os.path
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available tools for checking variants
VARIANT_CHECKERS = "peddy"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  variant_checking:
    tools_ngs_mapping: []  # optional, copied from ngs mapping config
    tools_variant_calling: []  # optional, copied from variant calling config
    path_variant_calling: ../variant_calling  # REQUIRED
    tools: ['peddy']
"""


class PeddyStepPart(BaseStepPart):
    """Compute variant statistics using peddy"""

    #: Step name
    name = "peddy"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{mapper}.{var_caller}.peddy.{index_ngs_library}/out/.done"
        self.log_path = (
            "work/{mapper}.{var_caller}.peddy.{index_ngs_library}/log/snakemake.filter.log"
        )

    @dictify
    def get_input_files(self, action):
        """Return path to pedigree input file"""
        # Validate action
        self._validate_action(action)
        yield "ped", os.path.realpath(
            "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        )
        tpl = (
            "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}"
        )
        key_ext = {"vcf": ".vcf.gz", "tbi": ".vcf.gz.tbi"}
        variant_calling = self.parent.sub_workflows["variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{mapper}.{var_caller}.peddy.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.peddy.{index_ngs_library}"
        )
        key_ext = {
            "background_pca": ".background_pca.json",
            "het_check": ".het_check.csv",
            "html": ".html",
            "ped_check": ".ped_check.csv",
            "ped": ".peddy.ped",
            "sex_check": ".sex_check.csv",
        }
        for key, ext in key_ext.items():
            yield key, prefix + ext

    def get_log_file(self, action):
        # Validate action
        self._validate_action(action)
        return self.log_path

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.

        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=1,
            time="10:00:00",  # 10 hours
            memory=f"{15 * 1024}M",
        )


class VariantCheckingWorkflow(BaseStep):
    """Perform germline variant checking"""

    #: Workflow name
    name = "variant_checking"

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

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
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((PeddyStepPart, WritePedigreeStepPart, LinkOutStepPart))
        # Register sub workflows
        self.register_sub_workflow("variant_calling", self.config["path_variant_calling"])
        # Copy over "tools" setting from ngs_mapping/variant_calling if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"]
        if not self.config["tools_variant_calling"]:
            self.config["tools_variant_calling"] = self.w_config["step_config"]["variant_calling"][
                "tools"
            ]

    @listify
    def get_result_files(self):
        """Return list of result files for the variant checking workflow"""
        yield from self._yield_peddy_results()

    def _yield_peddy_results(self):
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                for path in self.sub_steps["peddy"].get_output_files("run").values():
                    yield from expand(
                        path,
                        mapper=self.config["tools_ngs_mapping"]["dna"],
                        var_caller=self.config["tools_variant_calling"],
                        index_ngs_library=[pedigree.index.dna_ngs_library.name],
                    )
        if "gatk_hc_gvcf" in self.config["tools_variant_calling"]:
            for path in self.sub_steps["peddy"].get_output_files("run").values():
                yield from expand(
                    path,
                    mapper=self.config["tools_ngs_mapping"]["dna"],
                    var_caller=["gatk_hc_gvcf"],
                    index_ngs_library=["whole_cohort"],
                )

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_checking", "path_variant_calling"),
            "Path to variant calling not configured but required for variant checking",
        )

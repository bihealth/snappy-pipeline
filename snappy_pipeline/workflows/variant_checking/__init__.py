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

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow
from snappy_pipeline.workflows.variant_calling.model import ExpectedGermlineVariants

from .model import VariantChecking as VariantCheckingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available tools for checking variants
VARIANT_CHECKERS = "peddy"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = VariantCheckingConfigModel.default_config_yaml_string()


class PeddyStepPart(BaseStepPart):
    """Compute variant statistics using peddy"""

    #: Step name
    name = "peddy"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.variant_tool = str(self.parent.get_task_config("variant_calling").tool)
        self.base_path_out = "work/{index_ngs_library}/out/.done"
        self.log_path = "work/{index_ngs_library}/log/snakemake.filter.log"

    @dictify
    def get_input_files(self, action):
        """Return path to pedigree input file"""
        # Validate action
        self._validate_action(action)
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

        calling: ExpectedGermlineVariants = self.parent.get_upstream_paths(
            "variant_calling",
            library_name="{index_ngs_library}",
        )
        yield "vcf", calling.vcf
        yield "vcf_tbi", calling.vcf_tbi

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        # Validate action
        self._validate_action(action)
        prefix = "work/{index_ngs_library}/out/{index_ngs_library}"
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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
            runtime="10h",  # 10 hours
            mem=f"{15 * 1024}MB",
        )


class VariantCheckingWorkflow(BaseStep):
    """Perform germline variant checking"""

    #: Workflow name
    name = "variant_checking"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"germline"})): True}
    produces = [DataSignature(DataType.QC, frozenset({"pedigree_check"}))]
    config_model_class = VariantCheckingConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local variant-checking output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"ped_check": f"output/{lib}/out/{lib}.ped_check.csv"}

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
            previous_steps=(VariantCallingWorkflow, NgsMappingWorkflow),
            task_name=task_name,
            **kwargs,
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((PeddyStepPart, WritePedigreeStepPart, LinkOutStepPart))

        # Copy over "tools" setting from ngs_mapping/variant_calling if not set here

    @listify
    def get_result_files(self):
        """Return list of result files for the variant checking workflow"""
        yield from self._yield_peddy_results()

    def _yield_peddy_results(self):
        for index_library in self.output_entities:
            for path in self.sub_steps["peddy"].get_output_files("run").values():
                yield from expand(
                    path,
                    index_ngs_library=[index_library],
                )

# -*- coding: utf-8 -*-
"""Implementation of purity and ploidy checking for somatic NGS samples

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_purity_ploidy_estimate.rst

"""

import os
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand, touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.models import RelationshipDefinition

from .model import SomaticPurityPloidyEstimate as SomaticPurityPloidyEstimateConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Tools for estimating purity and ploidy.
PURITY_PLOIDY_TOOLS = "ascat"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticPurityPloidyEstimateConfigModel.default_config_yaml_string()

#: Extensions of output payload files
EXT_VALUES = ("_goodness_of_fit.txt", "_ploidy.txt", "_segments.txt", "_segments_raw.txt")


class AscatStepPart(BaseStepPart):
    """Estimation of purity and ploidy using ASCAT.

    Notes
    -----

    - Although we name the virtual probes "SNP${num}", they are not guaranteed
      to
    """

    # Name of the step.
    name = "ascat"

    # The actions for generating BAF/CNV files for the tumor and/nor normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "baf_tumor",
        "baf_normal",
        "cnv_tumor",
        "cnv_normal",
        "cnv_tumor_wes",
        "cnv_normal_wes",
        "run_ascat",
    )

    def __init__(self, parent):
        super().__init__(parent)

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def get_input_files(self, action):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))()

    def _get_input_files_baf_tumor(self):
        """Return input files for generating BAF file for the tumor."""

        def func(wildcards):
            ngs_mapping = self.parent.upstream("ngs_mapping")
            base_path = ("output/{tumor_library}/out/{tumor_library}").format(**wildcards)
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        return func

    def _get_input_files_baf_normal(self):
        """Return input files for generating BAF file for the normal."""

        def func(wildcards):
            ngs_mapping = self.parent.upstream("ngs_mapping")
            base_path = ("output/{normal_library}/out/{normal_library}").format(**wildcards)
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        return func

    def _get_input_files_cnv_tumor(self):
        """Return input files for generating BAF file for the tumor."""
        return self._get_input_files_baf_tumor()

    def _get_input_files_cnv_normal(self):
        """Return input files for generating CNV file for the normal."""
        return self._get_input_files_baf_normal()

    def _get_input_files_cnv_tumor_wes(self):
        """Return input files for generating CNV file from copywriter for tumor."""

        def func(wildcards):
            base_path = ("work/copywriter.{tumor_library}/out/copywriter.{tumor_library}").format(
                **wildcards
            )
            return {
                "bins": self.parent.upstream("somatic_targeted_seq_cnv_calling")(
                    base_path + "_bins.txt"
                )
            }

        return func

    def _get_input_files_cnv_normal_wes(self):
        """Return input files for generating CNV file from copywriter for normal."""

        def func(wildcards):
            df = self.parent.build_library_dataframe()
            normal_df = df[df["library_name"] == wildcards["normal_library"]]
            tumor_library = normal_df.iloc[0].get("library_name") if not normal_df.empty else None
            # Find tumor library that has this normal as matched_normal_lib
            if tumor_library is None:
                tumor_df = df[df["matched_normal_lib"] == wildcards["normal_library"]]
                if not tumor_df.empty:
                    tumor_library = tumor_df.iloc[0]["library_name"]
            base_path = ("work/copywriter.{tumor_library}/out/copywriter.{tumor_library}").format(
                tumor_library=tumor_library, **wildcards
            )
            return {
                "bins": self.parent.upstream("somatic_targeted_seq_cnv_calling")(
                    base_path + "_bins.txt"
                )
            }

        return func

    def _get_input_files_run_ascat(self):
        """Return input files for actually running ASCAT."""

        @dictify
        def func(wildcards):
            result = {
                "baf_tumor": (
                    "work/ascat_baf_tumor.{tumor_library}/out/ascat_baf_tumor.{tumor_library}.txt"
                ),
                "baf_normal": (
                    "work/ascat_baf_normal.{normal_library}/out/"
                    "ascat_baf_normal.{normal_library}.txt"
                ),
                "cnv_tumor": (
                    "work/ascat_cnv_tumor.{tumor_library}/out/ascat_cnv_tumor.{tumor_library}.txt"
                ),
                "cnv_normal": (
                    "work/ascat_cnv_normal.{normal_library}/out/"
                    "ascat_cnv_normal.{normal_library}.txt"
                ),
            }
            normal_library = self.get_normal_lib_name(wildcards)
            for key, value in result.items():
                yield key, value.format(normal_library=normal_library, **wildcards)

        return func

    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @staticmethod
    def _get_output_files_baf_tumor():
        """Return output files for generating BAF file for the tumor."""
        return {
            "txt": ("work/ascat_baf_tumor.{tumor_library}/out/ascat_baf_tumor.{tumor_library}.txt")
        }

    @staticmethod
    def _get_output_files_baf_normal():
        """Return output files for generating BAF file for the normal."""
        return {
            "txt": (
                "work/ascat_baf_normal.{normal_library}/out/ascat_baf_normal.{normal_library}.txt"
            )
        }

    @staticmethod
    def _get_output_files_cnv_tumor():
        """Return output files for generating BAF file for the tumor."""
        return {
            "txt": ("work/ascat_cnv_tumor.{tumor_library}/out/ascat_cnv_tumor.{tumor_library}.txt")
        }

    @staticmethod
    def _get_output_files_cnv_normal():
        """Return output files for generating CNV file for the normal."""
        return {
            "txt": (
                "work/ascat_cnv_normal.{normal_library}/out/ascat_cnv_normal.{normal_library}.txt"
            )
        }

    @dictify
    def _get_output_files_run_ascat(self):
        """Return output files for actually running ASCAT."""
        yield "done", touch("work/ascat.{tumor_library}/out/.done")
        infixes = ("goodness_of_fit", "ploidy", "segments", "segments_raw")
        for infix in infixes:
            path = ("work/ascat.{tumor_library}/out/{tumor_library}_%s.txt") % infix
            yield infix, path

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_baf_tumor(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "b_af_loci": self.config.ascat.b_af_loci,
            "reference_path": self.w_config.static_data_config.reference.path,
        }

    def _get_args_baf_normal(self, wildcards: Wildcards) -> dict[str, Any]:
        return self._get_args_baf_tumor(wildcards)

    def _get_args_cnv_tumor(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "b_af_loci": self.config.ascat.b_af_loci,
            "reference_path": self.w_config.static_data_config.reference.path,
            "tumor_library": wildcards.tumor_library,
        }

    def _get_args_cnv_normal(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "b_af_loci": self.config.ascat.b_af_loci,
            "reference_path": self.w_config.static_data_config.reference.path,
        }

    def _get_args_run_ascat(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"tumor_library": wildcards.tumor_library}

    def get_log_file(self, action):
        """Return path to log file"""
        # TODO: implement log option for actions `cnv_tumor_wes` and `cnv_normal_wes`.
        # Validate action
        self._validate_action(action)
        log_dict = {
            "baf_tumor": (
                "work/ascat_baf_tumor.{tumor_library}/log/ascat_baf_tumor.{tumor_library}.log"
            ),
            "baf_normal": (
                "work/ascat_baf_normal.{normal_library}/log/ascat_baf_normal.{normal_library}.log"
            ),
            "cnv_tumor": (
                "work/ascat_cnv_tumor.{tumor_library}/log/ascat_cnv_tumor.{tumor_library}.log"
            ),
            "cnv_normal": (
                "work/ascat_cnv_normal.{normal_library}/log/ascat_cnv_normal.{normal_library}.log"
            ),
            "run_ascat": ("work/ascat.{tumor_library}/log/ascat.{tumor_library}.log"),
        }
        return {"log": log_dict[action]}

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            runtime="2d",  # 2 days
            mem=f"{10 * 1024 * 8}MB",
        )


class SomaticPurityPloidyEstimateWorkflow(BaseStep):
    """Perform purity and ploidy estimation"""

    #: Workflow name
    name = "somatic_purity_ploidy_estimate"
    consumes = {DataSignature(DataType.VARIANTS, frozenset({"somatic", "cnv"})): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"purity_ploidy"}))]

    config_model_class = SomaticPurityPloidyEstimateConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local purity/ploidy output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"output/ascat.{lib}/out/.done"}

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
            previous_steps=(NgsMappingWorkflow,),
            task_name=task_name,
            **kwargs,
        )
        self.register_sub_step_classes((AscatStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the purity/ploidy estimation workflow."""
        tpl = os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}")
        for entity in self.output_entities:
            yield from expand(tpl, tumor_library=[entity], ext=EXT_VALUES)

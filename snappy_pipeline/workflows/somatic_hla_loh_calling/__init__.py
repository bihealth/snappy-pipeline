# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_hla_loh_calling`` step

This step allows for the detection of loss of heterzygosity for cancer samples
from whole genomes, exomes or large panels).
LOHHLA starts from the aligned reads, germline HLA calls and optionally purity
and plodiy estimates of a sample.

==========
Step Input
==========

``somatic_hla_loh_calling`` starts off the aligned reads, i.e. ``ngs_mapping``,
HLA calls from ``hla_calling`` and results from
``somatic_purity_ploidy_estimate``

===========
Step Output
===========

A report file.
"""

import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.hla_typing.model import ExpectedHlaTyping
from snappy_pipeline.models import RelationshipDefinition

from .model import SomaticHlaLohCalling as SomaticHlaLohCallingConfigModel

__author__ = "Clemens Messerschmidt <clemens.messerschmidt@bih-charite.de>"

#: Default configuration for the somatic_msi_calling step
DEFAULT_CONFIG = SomaticHlaLohCallingConfigModel.default_config_yaml_string()


class LohhlaStepPart(BaseStepPart):
    """Perform LOHHLA analysis of tumor/normal WES/WGS pairs."""

    #: Step name
    name = "lohhla"

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{tumor_library}/out/{tumor_library}{ext}"

    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)

        def input_function(wildcards):
            """Helper wrapper function"""
            normal_lib = self.get_normal_lib_name(wildcards)
            tumor_lib = wildcards.tumor_library
            normal: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=normal_lib
            )
            tumor: ExpectedAlignments = self.parent.get_upstream_paths(
                "ngs_mapping", library_name=tumor_lib
            )
            hla_typing: ExpectedHlaTyping = self.parent.get_upstream_paths(
                "hla_typing", library_name=normal_lib
            )
            return {
                "normal_bam": normal.bam,
                "normal_bai": normal.bai,
                "tumor_bam": tumor.bam,
                "tumor_bai": tumor.bai,
                "hla": hla_typing.txt,
            }

        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.tumor_library]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def get_output_files(self, action):
        """Return output files from LOHHLA"""
        # Validate action
        self._validate_action(action)
        return {"done": self.base_path_out.replace("{ext}", ".done")}

    @dictify
    def _get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = "work/{tumor_library}/log/{tumor_library}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"


class SomaticHlaLohCallingWorkflow(BaseStep):
    """Perform somatic hla loh calling"""

    #: Workflow name
    name = "somatic_hla_loh_calling"
    consumes = {DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})): True}
    produces = [DataSignature(DataType.TABULAR, frozenset({"hla_loh"}))]

    config_model_class = SomaticHlaLohCallingConfigModel

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
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        """Return local HLA LOH output paths for downstream consumers."""
        cls.require_signature(signature)
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"output/{lib}/out/{lib}.done"}

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
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((LohhlaStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        """Return list of result files for the somatic hla loh calling workflow."""
        log_exts = (
            ".log",
            ".log.md5",
            ".conda_info.txt",
            ".conda_info.txt.md5",
            ".conda_list.txt",
            ".conda_list.txt.md5",
        )
        for entity in self.output_entities:
            yield from expand(
                os.path.join("output", "{tumor_library}", "out", "{tumor_library}{ext}"),
                tumor_library=[entity],
                ext=".done",
            )
            yield from expand(
                os.path.join("output", "{tumor_library}", "log", "{tumor_library}{ext}"),
                tumor_library=[entity],
                ext=log_exts,
            )

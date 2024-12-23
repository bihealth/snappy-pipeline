# -*- coding: utf-8 -*-
"""Guess donor sex from NGS data (usually coverage of sex chromosomes vs autosomes)

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_guess_sex.rst

"""

from biomedsheets.shortcuts import GenericSampleSheet
from biomedsheets.io_tsv import EXTRACTION_TYPE_DNA
from snakemake.io import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.samplesheet import sample_sheets
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage

from .model import GuessSex as GuessSexModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = GuessSexModel.default_config_yaml_string()


class SamtoolsStepPart(BaseStepPart):
    """Estimation of donor sex using sample coverage over autosomes & sex chromosomes."""

    # Name of the step.
    name = "samtools"

    actions = ("run",)

    resource_usage = {"run": ResourceUsage(threads=2, time="02:00:00", memory=f"{4 * 1024 * 1}M")}

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg = getattr(self.config, self.name, {})

    def get_input_files(self, action: str):
        """Return input files"""
        # Validate action
        self._validate_action(action)

        def input_function(wildcards: Wildcards) -> dict[str, str]:
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(
                **wildcards
            )
            return {
                "bam": ngs_mapping(base_path + ".bam"),
                "bai": ngs_mapping(base_path + ".bam.bai"),
            }

        return input_function

    @dictify
    def get_output_files(self, action: str):
        # Validate action
        self._validate_action(action)
        tpl = "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}."
        output_files = {"table": tpl + "tsv", "decision": tpl + "txt"}
        for k, v in output_files.items():
            yield k, v
            yield k + "_md5", v + ".md5"

    def get_result_files(self) -> list[str]:
        tpl = "work/{mapper}.samtools.{library_name}/out/{mapper}.samtools.{library_name}."
        results = [tpl + "txt"]
        tpl = "work/{mapper}.samtools.{library_name}/log/{mapper}.samtools.{library_name}."
        for ext in ("conda_info.txt", "conda_list.txt", "log", "sh"):
            results.append(tpl + f"{ext}")
        return results

    @dictify
    def get_log_file(self, action: str):
        # Validate action
        self._validate_action(action)
        tpl = "work/{mapper}.samtools.{library_name}/log/{mapper}.samtools.{library_name}."
        for ext in ("conda_info.txt", "conda_list.txt", "log"):
            k = ext.replace(".txt", "")
            v = tpl + ext
            yield k, v
            yield k + "_md5", v + ".md5"
        k = "script"
        v = tpl + "sh"
        yield k, v
        yield k + "_md5", v + ".md5"

    def get_args(self, action: str) -> dict[str, str]:
        # Validate action
        self._validate_action(action)
        return dict(self.cfg)


class GuessSexWorkflow(BaseStep):
    """Guess donor sex from sample data"""

    #: Workflow name
    name = "guess_sex"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

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
            config_model_class=GuessSexModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        self.table = sample_sheets(self.sheets)
        self.table = self.table[self.table["extractionType"] == EXTRACTION_TYPE_DNA]
        self.register_sub_step_classes((SamtoolsStepPart, LinkOutStepPart))
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        for sub_step in self.sub_steps.values():
            tool = sub_step.name
            tool_config = getattr(self.config, tool, None)
            if tool_config is None:
                continue
            for row in self.table.itertuples():
                files_in_work = sub_step.get_result_files()
                for file_in_work in files_in_work:
                    file_in_output = file_in_work.format(
                        mapper=self.config.tool_ngs_mapping, library_name=row.Index
                    ).replace("work/", "output/", 1)
                    yield file_in_output
                    yield file_in_output + ".md5"

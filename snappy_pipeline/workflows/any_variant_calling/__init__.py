# -*- coding: utf-8 -*-
"""Abstract class for somatic & germline variant callers"""

import os

from snakemake.io import expand
from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import listify
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.common.samplesheet import filter_table_by_modality, sample_sheets

from .model import AnyVariantCalling as AnyVariantCallingConfigModel
from .model import VariantOrigin

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Extensions of log files (TODO: should go in a generic location)
LOG_EXT_VALUES = (
    ".log",
    ".log.md5",
    ".conda_info.txt",
    ".conda_info.txt.md5",
    ".conda_list.txt",
    ".conda_list.txt.md5",
)

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = AnyVariantCallingConfigModel.default_config_yaml_string()


class AnyVariantCallingWorkflow(BaseStep):
    """Placeholder class to facilitate generic variant annotation & filtration"""

    #: Workflow name
    name = "any_variant_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = GenericSampleSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        config_model_class=AnyVariantCallingConfigModel,
        previous_steps=(NgsMappingWorkflow,),
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=config_model_class,
            previous_steps=previous_steps,
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        self.table = filter_table_by_modality(sample_sheets(self.sheets), modality="dna")
        assert self.table.shape[1] > 0, "No valid samples"

        self.tools = self.config.tools
        assert self.tools, "No tool configured for variant calling"
        self.mapping_tools = set(self.w_config.step_config["ngs_mapping"].tools.dna)
        if self.config.tools_ngs_mapping:
            self.mapping_tools &= set(self.config.tools_ngs_mapping)
        assert self.mapping_tools, "No configured mapping tool"

    @listify
    def _get_result_only_files(self):
        """Returns the vcf files (vcf, vcf.tbi & checksums)"""
        name_pattern = "{mapper}.{caller}.{library}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.mapping_tools,
            caller=self.tools,
            ext=EXT_VALUES,
        )

    @listify
    def _get_log_only_files(self):
        """Returns the log files"""
        name_pattern = "{mapper}.{caller}.{library}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.mapping_tools,
            caller=self.tools,
            ext=LOG_EXT_VALUES,
        )

    def get_result_files(self):
        return self._get_result_only_files() + self._get_log_only_files()

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list."""
        match self.config.variant_origin:
            case VariantOrigin.SOMATIC:
                libraries = self.table[self.table["isTumor"]]
            case VariantOrigin.GERMLINE:
                libraries = self.table[~self.table["isTumor"]]
            case _:
                libraries = self.table
        for library in libraries["ngs_library"]:
            yield from expand(tpl, library=[library], **kwargs)

# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_neoepitope_prediction`` step

The somatic_neoepitope_prediction step allows for the prediction of neoepitopes from somatic
(small) variant calling results and a transcript database such as ENSEMBL.  Further, the step
allows for the binding prediction to a given set of HLA alleles.

.. note::

    Status: not implemented yet

==========
Step Input
==========

.. note:: TODO

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_neoepitope_prediction.rst

"""

from collections import OrderedDict
import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_neoepitope_prediction:
    path_somatic_variant_annotation: ../somatic_variant_annotation  # REQUIRED
    path_rna_ngs_mapping: ../ngs_mapping
    tools_somatic_variant_annotation: ["vep"]
    tools_rna_mapping: [] # deafult to those configured for ngs_mapping
    tools_ngs_mapping: [] # deafult to those configured for ngs_mapping
    tools_somatic_variant_calling: [] # deafult to those configured for somatic_variant_calling
    preparation:
      format: 'snappy_custom' # REQUIRED - The file format of the expression file to process. (stringtie, kallisto, cufflinks, snappy_custom, custom)
      # Use `custom` to process file formats not explicitly supported.
      # The `custom` option requires the use of the --id-column and --expression-column arguments.
      expression_file: '' # expression file path. If the format is not snappy_custom
      path_features: '' # Gencode file path, required for star and snappy format
      mode: 'gene'  # REQUIRED - Determine whether the expression file contains gene or transcript TPM values.
      full_vep_annotation: True # The somatic_variant_annotation has been done in fully annotation mode
      id_column:  # Gene/Transcript id column name. Need when using the `custom` format.
      expression_column:  # Expression column name. Need when using the `custom` format.
      ignore_ensembl_id_version: True #Ignore the ensemble id version
      max_depth: 4000
"""


class NeoepitopePreparationStepPart(BaseStepPart):
    """
    Preparation VCF file for pvactool
    """

    #: Step name
    name = "neoepitope_preparation"
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @dictify
    def get_input_files(self, action):
        """Return path to somatic variant annotated file"""
        # It only works with vep now.
        # Validate action
        self._validate_action(action)
        tpl = (
            "output/{mapper}.{var_caller}.{anno_caller}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{anno_caller}.{tumor_library}"
        )
        if self.config["preparation"]["full_vep_annotation"]:
            key_ext = {"vcf": ".full.vcf.gz", "vcf_tbi": ".full.vcf.gz.tbi"}
        else:
            key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        for key, ext in key_ext.items():
            yield key, variant_annotation(tpl + ext)
        # Getting appropriate rna library
        if self.config["preparation"]["format"] == "snappy_custom":
            for sheet in filter(is_not_background, self.parent.sheets):
                for donor in sheet.bio_entities.values():
                    for biosample in donor.bio_samples.values():
                        if biosample.extra_infos["isTumor"]:
                            for test_sample in biosample.test_samples.values():
                                if test_sample.extra_infos["extractionType"] == "RNA":
                                    for lib in test_sample.ngs_libraries.values():
                                        rna_library = lib.name
                                        rna_tpl = (
                                            "output/{mapper}.{library_name}/out/{mapper}.{library_name}"
                                        ).format(
                                            mapper="star",
                                            library_name=rna_library,
                                        )
                                        ext = {"expression", "bam", "bai"}
                                        yield "expression", ngs_mapping(rna_tpl + ".GeneCounts.tab")
                                        yield "bam", ngs_mapping(rna_tpl + ".bam")
                                        yield "bai", ngs_mapping(rna_tpl + ".bam.bai")

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        if self.config["preparation"]["mode"] == "gene":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
            )
        elif self.config["preparation"]["mode"] == "transcript":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}"
            )
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        """Return mapping of log files."""
        # Validate action
        self._validate_action(action)
        if self.config["preparation"]["mode"] == "gene":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/log/"
                "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
            )
        elif self.config["preparation"]["mode"] == "transcript":
            prefix = (
                "work/{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}/log/"
                "{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}"
            )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_resource_usage(self, action):
        self._validate_action(action)
        mem_mb = 6 * 1024  # 6GB
        return ResourceUsage(
            threads=2,
            time="1:00:00",  # 1 hour
            memory=f"{mem_mb}M",
        )

    def get_params(self, action):
        self._validate_action(action)
        return {
            "max_depth": self.config["preparation"]["max_depth"],
            "format": self.config["preparation"]["format"],
            "expression_column": self.config["preparation"]["expression_column"],
            "id_column": self.config["preparation"]["id_column"],
            "ignore_ensembl_id_version": self.config["preparation"]["ignore_ensembl_id_version"],
            "mode": self.config["preparation"]["mode"],
        }


class SomaticNeoepitopePredictionWorkflow(BaseStep):
    """Perform neoepitope prediction workflow"""

    name = "somatic_neoepitope_prediction"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
            (
                SomaticVariantCallingWorkflow,
                SomaticVariantAnnotationWorkflow,
                NgsMappingWorkflow,
            ),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((NeoepitopePreparationStepPart, LinkOutStepPart))
        self.register_sub_workflow(
            "somatic_variant_annotation",
            self.config["path_somatic_variant_annotation"],
        )
        self.register_sub_workflow(
            "ngs_mapping",
            self.config["path_rna_ngs_mapping"],
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.config["tools_ngs_mapping"]:
            self.config["tools_ngs_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "dna"
            ]
        if not self.config["tools_somatic_variant_calling"]:
            self.config["tools_somatic_variant_calling"] = self.w_config["step_config"][
                "somatic_variant_calling"
            ]["tools"]
        if not self.config["tools_somatic_variant_annotation"]:
            self.config["tools_somatic_variant_annotation"] = ["vep"]
        if not self.config["tools_rna_mapping"]:
            self.config["tools_rna_mapping"] = self.w_config["step_config"]["ngs_mapping"]["tools"][
                "rna"
            ]
        if not self.config["preparation"]["path_features"]:
            self.config["preparation"]["path_features"] = self.w_config["static_data_config"][
                "features"
            ]["path"]

    @listify
    def get_result_files(self):
        callers = set(self.config["tools_somatic_variant_calling"])
        anno_callers = set(self.config["tools_somatic_variant_annotation"])
        if self.config["preparation"]["mode"] == "gene":
            name_pattern = "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library.name}"
        elif self.config["preparation"]["mode"] == "transcript":
            name_pattern = "{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            anno_caller=anno_callers,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.config["tools_ngs_mapping"],
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            anno_caller=anno_callers,
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                    or not sample_pair.tumor_sample.rna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl,
                    tumor_library=[sample_pair.tumor_sample.dna_ngs_library],
                    **kwargs,
                )

    def check_config(self):
        """Check that path_somatic_variant_annotation, preparation/mode and preparation/format are present in the configuration"""
        self.ensure_w_config(
            (
                "step_config",
                "somatic_neoepitope_prediction",
                "path_somatic_variant_annotation",
            ),
            "Path to variant (directory of vcf files) not configured but required for somatic neoepitope prediction",
        )

        self.ensure_w_config(
            ("step_config", "somatic_neoepitope_prediction", "preparation", "mode"),
            "The mode is required for adding gene expression data to somatic variant annotated file",
        )

        self.ensure_w_config(
            ("step_config", "somatic_neoepitope_prediction", "preparation", "format"),
            "Format is required to determine which tool generated the TPM values. ",
        )

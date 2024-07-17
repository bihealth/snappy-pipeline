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
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from .model import SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel


__author__ = "Pham Gia Cuong"
__email__ = "pham.gia-cuong@bih-charite.de"

#: Extensions of files to create as main payload
PREPARE_EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
PREDICT_EXT_VALUES = (
    ".all_epitopes.tsv",
    ".filtered.tsv",
    ".all_epitopes.tsv.md5",
    ".filtered.tsv.md5",
)
#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticNeoepitopePredictionConfigModel.default_config_yaml_string()


class PvacSeqStepPart(BaseStepPart):
    """
    Preparation VCF file for pvactool
    """

    #: Step name
    name = "pvacseq"

    #: Actions
    actions = ("install", "prepare", "predict")

    #: Resources
    resource_usage = {
        "install": ResourceUsage(
            threads=2,
            time="03:00:00",
            memory="6G",
        ),
        "prepare": ResourceUsage(
            threads=1,
            time="01:00:00",
            memory="6G",
        ),
        "predict": ResourceUsage(
            threads=4,
            time="4:00:00",
            memory="30G",
        ),
    }

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

    def get_input_files(self, action):
        self._validate_action(action)
        # if action == "install":
        #     return {"container": "work/containers/out/pvactools.simg"}
        if action == "prepare":
            return self._get_input_files_prepare
        if action == "predict":
            return self._get_input_files_predict

    @dictify
    def _get_input_files_predict(self, wildcards):
        yield "container", "work/containers/out/pvactools.simg"
        prepare_tpl = (
            "work/prepare/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}.vcf.gz"
        )
        yield "combine_vcf", prepare_tpl
        hla_typing = self.parent.sub_workflows["hla_typing"]
        hla_tpl = "output/optitype.{ngs_library}/out/optitype.{ngs_library}.txt"
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if sample_pair.tumor_sample.dna_ngs_library.name == wildcards.tumor_library:
                    yield (
                        "hla_normal_dna",
                        hla_typing(
                            (hla_tpl).format(
                                ngs_library=sample_pair.normal_sample.dna_ngs_library.name,
                            )
                        ),
                    )
                    yield (
                        "hla_tumor_dna",
                        hla_typing(
                            (hla_tpl).format(
                                ngs_library=sample_pair.tumor_sample.dna_ngs_library.name,
                            )
                        ),
                    )
                    yield (
                        "hla_tumor_rna",
                        hla_typing(
                            (hla_tpl).format(
                                ngs_library=sample_pair.tumor_sample.rna_ngs_library.name,
                            )
                        ),
                    )

    @dictify
    def _get_input_files_prepare(self, wildcards):
        """Return path to somatic variant annotated file"""
        # It only works with vep now.
        # Validate action
        # self._validate_action(action)
        # yield "container","work/containers/out/pvactools.simg"
        tpl = (
            "output/{mapper}.{var_caller}.{anno_caller}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{anno_caller}.{tumor_library}"
        )
        if self.config["preparation"]["full_vep_annotation"]:
            vep_key_ext = {"vcf": ".full.vcf.gz", "vcf_tbi": ".full.vcf.gz.tbi"}
        else:
            vep_key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        if self.config["preparation"]["format"] == "snappy_custom":
            for key, ext in vep_key_ext.items():
                yield key, variant_annotation(tpl + ext)

        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        star_key_ext = {"expression": ".GeneCounts.tab", "bam": ".bam"}
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if sample_pair.tumor_sample.dna_ngs_library.name == wildcards.tumor_library:
                    rna_tpl = "output/{mapper}.{rna_library}/out/{mapper}.{rna_library}"
                    for key, ext in star_key_ext.items():
                        yield (
                            key,
                            ngs_mapping(
                                (rna_tpl + ext).format(
                                    mapper=self.w_config.step_config["ngs_mapping"].tools.rna[0],
                                    rna_library=sample_pair.tumor_sample.rna_ngs_library.name,
                                )
                            ),
                        )

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        # Validate action
        self._validate_action(action)
        if action == "install":
            yield "container", "work/containers/out/pvactools.simg"
        if action == "prepare":
            path_prepare = self._get_output_files_prepare()
            yield from path_prepare.items()
        if action == "predict":
            path_predict = self._get_output_files_predict()
            yield from path_predict.items()

    @dictify
    def _get_output_files_prepare(self):
        if self.config["preparation"]["mode"] == "gene":
            prefix = (
                "work/prepare/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
            )
        elif self.config["preparation"]["mode"] == "transcript":
            prefix = (
                "work/prepare/{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}/out/"
                "{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}"
            )
        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_output_files_predict(self):
        key_ext = {"all_epitopes": ".all_epitopes.tsv", "filtered_epitopes": ".filtered.tsv"}
        prefix = "work/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}.epitopes/out/MHC_Class_I/{tumor_library}"
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        """Return mapping of log files."""
        # Validate action
        self._validate_action(action)
        if action == "install":
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
            )
            for key, ext in key_ext:
                yield key, "work/containers/log/pvactool_install" + ext

        if action == "prepare":
            if self.config["preparation"]["mode"] == "gene":
                prefix = (
                    "work/prepare/{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}/log/"
                    "{mapper}.{var_caller}.{anno_caller}.GX.{tumor_library}"
                )
            elif self.config["preparation"]["mode"] == "transcript":
                prefix = (
                    "work/preprare/{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}/log/"
                    "{mapper}.{var_caller}.{anno_caller}.TX.{tumor_library}"
                )
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
            )
            for key, ext in key_ext:
                yield key, prefix + ext

        if action == "predict":
            prefix = (
                "work/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library}.epitopes/"
                "log/{tumor_library}"
            )
            key_ext = (
                ("log", ".log"),
                ("conda_info", ".conda_info.txt"),
                ("conda_list", ".conda_list.txt"),
            )
            for key, ext in key_ext:
                yield key, prefix + ext

    def get_params(self, action):
        self._validate_action(action)
        if action == "prepare":
            return getattr(self, "_get_params_run_prepare")
        if action == "predict":
            return getattr(self, "_get_params_run_predict")

    def _get_params_run_prepare(self, wildcards):
        d = {
            "max_depth": self.config["preparation"]["max_depth"],
            "format": self.config["preparation"]["format"],
            "expression_column": self.config["preparation"]["expression_column"],
            "id_column": self.config["preparation"]["id_column"],
            "ignore_ensembl_id_version": self.config["preparation"]["ignore_ensembl_id_version"],
            "mode": self.config["preparation"]["mode"],
        }
        return d

    def _get_params_run_predict(self, wildcards):
        d = {
            "epitope_length": ",".join(
                map(str, self.config["prediction"]["CLASS_I_EPITOPE_LENGTH"])
            ),
            "algorithms": " ".join(self.config["prediction"]["algorithms"]),
            "BINDING_THRESHOLD": self.config["prediction"]["BINDING_THRESHOLD"],
            "percentile_threshold": self.config["prediction"]["percentile_threshold"],
            "allele_specific_binding_thresholds": self.config["prediction"][
                "allele_specific_binding_thresholds"
            ],
            "aggregate_inclusion_binding_threshold": self.config["prediction"][
                "aggregate_inclusion_binding_threshold"
            ],
            "netmhc_stab": self.config["prediction"]["netmhc_stab"],
            "NET_CHOP_THRESHOLD": self.config["prediction"]["NET_CHOP_THRESHOLD"],
            "PROBLEMATIC_AMINO_ACIDS": self.config["prediction"]["PROBLEMATIC_AMINO_ACIDS"],
            # 'run_reference_proteome_similarity': self.config["prediction"]["run_reference_proteome_similarity"],
            "FAST_SIZE": self.config["prediction"]["FAST_SIZE"],
            "exclude_NAs": self.config["prediction"]["exclude_NAs"],
            "NORMAL_COV": self.config["prediction"]["NORMAL_COV"],
            "TDNA_COV": self.config["prediction"]["TDNA_COV"],
            "TRNA_COV": self.config["prediction"]["TRNA_COV"],
            "NORMAL_VAF": self.config["prediction"]["NORMAL_VAF"],
            "maximum_transcript_support_level": self.config["prediction"][
                "maximum_transcript_support_level"
            ],
            "pass_only": self.config["prediction"]["pass_only"],
            "tumor_purity": self.config["prediction"]["tumor_purity"],
        }
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if sample_pair.tumor_sample.dna_ngs_library.name == wildcards.tumor_library:
                    d["normal_sample"] = sample_pair.normal_sample.dna_ngs_library.name
        return d


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
            config_model_class=SomaticNeoepitopePredictionConfigModel,
            previous_steps=(
                SomaticVariantCallingWorkflow,
                SomaticVariantAnnotationWorkflow,
                NgsMappingWorkflow,
                HlaTypingWorkflow,
            ),
        )
        # Register sub step classes so the sub steps are available
        config = self.config
        self.register_sub_step_classes((PvacSeqStepPart, LinkOutStepPart))
        self.register_sub_workflow(
            "somatic_variant_annotation",
            config["path_somatic_variant_annotation"],
        )
        self.register_sub_workflow(
            "ngs_mapping",
            config["path_rna_ngs_mapping"],
        )
        self.register_sub_workflow(
            "hla_typing",
            config["path_hla_typing"],
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not config.tools_ngs_mapping:
            config.tools_ngs_mapping = self.w_config["step_config"]["ngs_mapping"].tools.dna
        if not config.tools_somatic_variant_calling:
            config.tools_somatic_variant_calling = self.w_config["step_config"][
                "somatic_variant_calling"
            ].tools
        if not config.tools_somatic_variant_annotation:
            config.tools_somatic_variant_annotation = ["vep"]
        if not config.tools_rna_mapping:
            config.tools_rna_mapping = self.w_config["step_config"]["ngs_mapping"].tools.rna
        if not config.preparation.path_features:
            config.preparation.path_features = self.w_config["static_data_config"].features.path

    @listify
    def get_result_files(self):
        cfg_mode = self.config["preparation"]["mode"]
        callers = set(self.config["tools_somatic_variant_calling"])
        anno_callers = set(self.config["tools_somatic_variant_annotation"])
        predict_tpl = "output/predict/{mapper}.{var_caller}.{anno_caller}.{mode}.{tumor_library.name}.epitopes/"
        yield from self._yield_result_files_matched(
            predict_tpl + "out/MHC_Class_I/{tumor_library.name}" + "{ext}",
            mapper=self.config["tools_ngs_mapping"],
            var_caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            anno_caller=anno_callers,
            mode="GX" if cfg_mode == "gene" else "TX",
            ext=PREDICT_EXT_VALUES,
        )

        yield from self._yield_result_files_matched(
            predict_tpl + "log/{tumor_library.name}" + "{ext}",
            mode="GX" if cfg_mode == "gene" else "TX",
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

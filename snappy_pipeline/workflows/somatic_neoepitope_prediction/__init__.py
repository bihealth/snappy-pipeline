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

import pandas as pd

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand, Wildcards, InputFiles

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.models.common import ExtractionType
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.common.samplesheet import sample_sheets, tumor_to_normal_mapping
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.gene_expression_quantification import (
    GeneExpressionQuantificationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from .model import SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel
from .model import PVACseq as PVACseqModel


__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

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
    actions = ("install", "pileup", "combine", "predict")

    #: Resources
    resource_usage = {
        "install": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "pileup": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "combine": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "predict": ResourceUsage(
            threads=4,
            time="23:59:59",
            memory="64G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: PVACseqModel = self.config.pVACseq
        self.annotated_template = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            self.annotated_template += ".filtered"
        self.annotated_template += ".{tumor_dna}"

    def get_input_files(self, action):
        match action:
            case "install":
                return lambda x: {}
            case "pileup":
                return self._get_input_files_pileup
            case "combine":
                return self._get_input_file_combine
            case "predict":
                return self._get_input_files_predict
            case _:
                self._validate_action(action)

    def _get_input_files_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {}

        tpl = "output/{{mapper}}.{library}/out/{{mapper}}.{library}.bam".format(
            library=self.parent.tumor_rna[wildcards.tumor_dna]
        )
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        input_files["bam"] = ngs_mapping(tpl)

        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.annotated_template)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        input_files["loci"] = annotation(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path

        return input_files

    def _get_input_file_combine(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {}

        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.annotated_template)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        input_files["annotated"] = annotation(tpl)

        if self.config.pileup.enabled and wildcards.tumor_dna in self.parent.tumor_rna:
            tpl = "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.annotated_template)
            input_files["pileup"] = tpl

        if self.config.quantification.enabled and (
            library := self.parent.tumor_rna.get(wildcards.tumor_dna, None)
        ):
            tpl = f"output/{{mapper}}.{library}/out/{{mapper}}.{library}.gene.sf"
            quantification = self.parent.sub_workflows["gene_expression_quantification"]
            input_files["gene_tpms"] = quantification(tpl)

            tpl = f"output/{{mapper}}.{library}/out/{{mapper}}.{library}.transcript.sf"
            input_files["transcript_tpms"] = quantification(tpl)

        return input_files

    @dictify
    def _get_input_files_predict(self, wildcards: Wildcards):
        if self.cfg.path_container:
            yield "container", self.cfg.path_container
        else:
            yield "container", "work/containers/out/pvactools.sif"

        yield "vcf", "work/{tpl}/out/{tpl}.combined.vcf.gz".format(tpl=self.annotated_template)

        hla_typing = self.parent.sub_workflows["hla_typing"]
        library = wildcards.tumor_dna
        tpl = f"{{typer}}.{library}"
        yield "hla_tumor_dna", hla_typing("output/" + tpl + "/out/" + tpl + ".txt")
        if library := self.parent.tumor_dna.get(wildcards.tumor_dna, None):
            tpl = f"{{typer}}.{library}"
            yield "hla_normal_dna", hla_typing("output/" + tpl + "/out/" + tpl + ".txt")
        if library := self.parent.tumor_rna.get(wildcards.tumor_dna, None):
            tpl = f"{{typer}}.{library}"
            yield "hla_tumor_rna", hla_typing("output/" + tpl + "/out/" + tpl + ".txt")

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.cfg.run_reference_proteome_similarity and self.cfg.peptide_fasta:
            yield "peptides", self.cfg.peptide_fasta

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        match action:
            case "install":
                return {"container": "work/containers/out/pvactools.sif"}
            case "pileup":
                return {
                    "vcf": "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.annotated_template)
                }
            case "combine":
                return {
                    "vcf": "work/{tpl}/out/{tpl}.combined.vcf.gz".format(
                        tpl=self.annotated_template
                    )
                }
            case "predict":
                return {"done": "work/{tpl}/neoepitopes/.done".format(tpl=self.annotated_template)}
            case _:
                self._validate_action(action)

    @dictify
    def get_log_file(self, action):
        """Return mapping of log files."""
        self._validate_action(action)
        if action == "install":
            tpl = "work/containers/log/pvactools"
        else:
            tpl = "work/{tpl}/log/{action}".format(tpl=self.annotated_template, action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    def get_args(self, action):
        match action:
            case "pileup":
                return self._get_args_pileup
            case "combine":
                return self._get_args_combine
            case "predict":
                return self._get_args_predict
            case _:
                self._validate_action(action)

    def _get_sample_names(self, wildcards: Wildcards) -> dict[str, str]:
        args = {}
        tumor = self.parent.sample_table[
            self.parent.sample_table["ngs_library"] == wildcards.tumor_dna
        ]
        args["tumor_sample"] = tumor["bio_sample"].iat[0]
        args["tumor_library"] = wildcards.tumor_dna
        if normal_dna := self.parent.tumor_dna.get(wildcards.tumor_dna, None):
            normal = self.parent.sample_table[self.parent.sample_table["ngs_library"] == normal_dna]
            args["normal_sample"] = normal["bio_sample"].iat[0]
            args["normal_library"] = normal_dna
        return args

    @staticmethod
    def _extra_args_flags(args: dict[str, Any]):
        for k in list(args.keys()):
            v = args[k]
            if isinstance(v, bool):
                args.pop(k)
                if v:
                    yield f"--{k.replace('_', '-')}"

    @staticmethod
    def _extra_args_lists(args: dict[str, Any], sep=",") -> dict[str, Any]:
        for k in list(args.keys()):
            v = args[k]
            if isinstance(v, list):
                if v:
                    args[k] = sep.join(map(str, v))
                else:
                    del args[k]
        return args

    @staticmethod
    def _group_extra_args(args: dict[str, Any]):
        for k, v in args.items():
            k = k.replace("_", "-")
            if isinstance(v, str):
                v = "'" + v + "'"
            yield f"--{k} {v}"

    def _get_args_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.pileup.model_dump(by_alias=True))

        del args["enabled"]
        del args["path_ngs_mapping"]
        del args["tool_rna_mapping"]
        if args["baq"] is None:
            del args["baq"]

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))

        return {
            "tumor_sample": self._get_sample_names(wildcards)["tumor_sample"],
            "extra_args": extra_args.strip(),
        }

    def _get_args_combine(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.quantification.model_dump(by_alias=True))

        del args["enabled"]
        del args["path_gene_expression_quantification"]
        args["format"] = args.pop("tool_gene_expression_quantification")

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))
        return {"extra_args": extra_args.strip()} | self._get_sample_names(wildcards)

    def _get_args_predict(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["peptide_fasta"]
        del args["genes_of_interest_file"]

        if isinstance(args["algorithms"], list):
            args["algorithms"] = ",".join(args["algorithms"])

        if args["percentile_threshold"] is None:
            del args["percentile_threshold"]
            del args["percentile_threshold_strategy"]
        else:
            del args["binding_threshold"]

        if args["maximum_transcript_support_level"] is None:
            del args["maximum_transcript_support_level"]

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))

        hla_types = set([])
        for i in ("hla_tumor_dna", "hla_normal_dna", "hla_tumor_rna"):
            if fn := getattr(input, i, None):
                with open(fn, "rt") as f:
                    for line in f:
                        hla_types += "HLA-" + line.strip()

        samples = self._get_sample_names(wildcards)

        return {
            "normal_sample": samples["normal_sample"],
            "tumor_sample": samples["tumor_sample"],
            "alleles": sorted(hla_types),
            "extra_args": extra_args.strip(),
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
        previous_steps: list[BaseStep] = [SomaticVariantAnnotationWorkflow, HlaTypingWorkflow]
        if config["step_config"]["somatic_neoepitope_prediction"]["pileup"]["enabled"]:
            previous_steps.append(NgsMappingWorkflow)
        if config["step_config"]["somatic_neoepitope_prediction"]["quantification"]["enabled"]:
            previous_steps.append(GeneExpressionQuantificationWorkflow)

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticNeoepitopePredictionConfigModel,
            previous_steps=previous_steps,
        )

        self.register_sub_step_classes((PvacSeqStepPart, LinkOutStepPart))

        self.register_sub_workflow(
            "somatic_variant_annotation",
            self.config.path_somatic_variant_annotation,
        )
        self.register_sub_workflow(
            "hla_typing",
            self.config.path_hla_typing,
        )
        if self.config.pileup.enabled:
            self.register_sub_workflow(
                "ngs_mapping",
                self.config.pileup.path_ngs_mapping,
            )
        if self.config.quantification.enabled:
            self.register_sub_workflow(
                "gene_expression_quantification",
                self.config.quantification.path_gene_expression_quantification,
            )

        self.sample_table: pd.DataFrame = sample_sheets(self.sheets)
        assert "extractionType" in self.sample_table.columns, (
            "'extractionType' missing from sample sheet"
        )
        self.tumor_dna = tumor_to_normal_mapping(
            self.sample_table[self.sample_table["extractionType"] == ExtractionType.DNA]
        )
        self.tumor_rna = self._dna_to_rna_mapping(self.sample_table)
        if self.config.pileup.enabled or self.config.quantification.enabled:
            assert any(map(lambda lib: lib in self.tumor_rna.keys(), self.tumor_dna.keys())), (
                "No tumor sample with somatic variant has expression data"
            )

    def _dna_to_rna_mapping(self, sample_table: pd.DataFrame) -> dict[str, str]:
        dna = sample_table[sample_table["extractionType"] == ExtractionType.DNA]
        rna = sample_table[sample_table["extractionType"] == ExtractionType.RNA]
        dna_rna_map = dna[["ngs_library", "bio_entity", "bio_sample"]].merge(
            rna[["ngs_library", "bio_entity", "bio_sample"]], on=["bio_entity", "bio_sample"]
        )
        return pd.Series(
            dna_rna_map.ngs_library_y.values, index=dna_rna_map.ngs_library_x.values
        ).to_dict()

    @listify
    def get_result_files(self):
        mappers = self.w_config.step_config["ngs_mapping"]["tools"]["dna"]
        callers = self.w_config.step_config["somatic_variant_calling"]["tools"]
        annotators = self.w_config.step_config["somatic_variant_annotation"]["tools"]
        tpl = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            tpl += ".filtered"
        tpl += ".{tumor_dna}"
        for tumor_dna in self.tumor_dna.keys():
            yield from expand(
                "output/" + tpl + "/neoepitopes/.done",
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                tumor_dna=[tumor_dna],
            )
            yield from expand(
                "output/" + tpl + "/log/predict.log",
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                tumor_dna=[tumor_dna],
            )

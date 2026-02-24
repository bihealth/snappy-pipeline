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

import json
import os
import re

import pandas as pd

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand, Wildcards, InputFiles

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.models.common import ExtractionType
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.common.samplesheet import sample_sheets, tumor_to_normal_mapping
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from .model import SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel
from .model import PVACseq as PVACseqModel
from .model import PVACfuse as PVACfuseModel
from .model import PVACsplice as PVACspliceModel
from .model import NetChopMethod


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

#: MHC class I & II genes & pseudo-genes
CLASS_I = ("A", "B", "C")  # , "E", "F", "G", "H", "J", "K", "L", "N", "S", "T", "U", "V", "W", "Y")
CLASS_II = ("DPA1", "DPB1", "DPA2", "DPB2", "DQA1", "DQB1")
# , "DRB1", "DRB2", "DRB3", "DRB4", "DRB5", "DRB7")

EXTRACTION_TYPES = ("dna", "rna")
MHC_CLASSES = ("class_i", "class_ii")


class UnsupportedProtocolStrand(Exception):
    pass


class PvacToolsStepPart(BaseStepPart):
    """Generic stuff for pVACtools modules"""

    #: Step name
    name = "pvactools"

    actions = ("install",)
    require_rna: bool = False

    CLASS_I_PATTERN = re.compile(
        r"^(?P<valid>({})\*[0-9]+:[0-9]+N?)(?P<suppl>.*)$".format("|".join(CLASS_I))
    )
    CLASS_II_PATTERN = re.compile(
        r"^(?P<valid>({})\*[0-9]+:[0-9]+N?)(?P<suppl>.*)$".format("|".join(CLASS_II))
    )

    def __init__(self, parent):
        super().__init__(parent)
        self.annotated_template = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            self.annotated_template += ".filtered"
        self.annotated_template += ".{tumor_dna}"

        self.hla_tools = {}
        for extraction_type in EXTRACTION_TYPES:
            for mhc_class in MHC_CLASSES:
                if tool := self.config.tools_hla_typing.get(extraction_type, {}).get(
                    mhc_class, None
                ):
                    if extraction_type not in self.hla_tools:
                        self.hla_tools[extraction_type] = {}
                    if (
                        mapper := self.w_config.step_config.get("hla_typing")
                        .get(tool)
                        .get("mapper", None)
                    ):
                        self.hla_tools[extraction_type][mhc_class] = f"{mapper}.{tool}"
                    else:
                        self.hla_tools[extraction_type][mhc_class] = tool

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def get_output_files(self, action):
        if action == "install":
            return {"container": "work/containers/out/pvactools.sif"}
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def _get_output_files_run(self):
        return {
            "done": "work/{tpl}/out/{name}.done".format(
                tpl=self.annotated_template, name=self.name
            ),
            "path": "work/{tpl}/{name}/.done".format(tpl=self.annotated_template, name=self.name),
        }

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def get_log_file(self, action):
        """Return mapping of log files."""
        if action == "install":
            return f"work/containers/log/{self.name}.log"
        self._validate_action(action)
        tpl = "work/{tpl}/log/{action}".format(tpl=self.annotated_template, action=action)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = tpl + ext
            log_files[key + "_md5"] = log_files[key] + ".md5"
        return log_files

    @listify
    def get_result_files(self):
        mappers = self.w_config.step_config["ngs_mapping"]["tools"]["dna"]
        callers = self.w_config.step_config["somatic_variant_calling"]["tools"]
        annotators = self.w_config.step_config["somatic_variant_annotation"]["tools"]
        tpl = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            tpl += ".filtered"
        tpl += ".{tumor_dna}"
        for tumor_dna in self.parent.tumor_dna.keys():
            if self.require_rna and self.parent.tumor_rna.get(tumor_dna, None) is None:
                continue
            yield from expand(
                "output/" + tpl + f"/out/{self.name}.done",
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                tumor_dna=[tumor_dna],
            )
            yield from expand(
                "output/" + tpl + f"/log/{self.name}.log",
                mapper=mappers,
                caller=callers,
                annotator=annotators,
                tumor_dna=[tumor_dna],
            )

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

    @listify
    def _get_hla_files(self, wildcards: Wildcards):
        hla_typing = self.parent.sub_workflows["hla_typing"]
        tumor_dna = wildcards.tumor_dna
        normal_dna = self.parent.tumor_dna.get(tumor_dna, None)
        tumor_rna = self.parent.tumor_rna.get(tumor_dna, None)

        input_files = []
        for mhc_class in MHC_CLASSES:
            if tool := self.hla_tools.get("dna", {}).get(mhc_class, None):
                tpl = "output/{tool}.{library_name}/out/{tool}.{library_name}.json"
                input_files.append(tpl.format(tool=tool, library_name=tumor_dna))
                if normal_dna:
                    input_files.append(tpl.format(tool=tool, library_name=normal_dna))
        if tumor_rna:
            for mhc_class in MHC_CLASSES:
                if tool := self.hla_tools.get("rna", {}).get(mhc_class, None):
                    tpl = "output/{tool}.{library_name}/out/{tool}.{library_name}.json"
                    input_files.append(tpl.format(tool=tool, library_name=tumor_rna))

        for f in input_files:
            yield hla_typing(f)

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

    @classmethod
    def _read_hla_values(cls, hla_typing_files: list[str], mhc_class: str) -> list[str]:
        if mhc_class == MHC_CLASSES[0]:
            loci = CLASS_I
            pattern = cls.CLASS_I_PATTERN
            prefix = "HLA-"
        elif mhc_class == MHC_CLASSES[1]:
            loci = CLASS_II
            pattern = cls.CLASS_II_PATTERN
            prefix = ""
        else:
            raise MissingConfiguration(f"Unknown MHC class {mhc_class}")

        hla_types = []
        for fn in hla_typing_files:
            with open(fn, "rt") as f:
                calls: dict[str, Any] = json.load(f)
            for locus, alleles in calls.items():
                if locus in loci:
                    for allele in alleles:
                        m = pattern.match(allele)
                        if m:
                            hla_types.append(prefix + m.group("valid"))

        return list(set(hla_types))

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        if action == "install":
            return ResourceUsage(threads=1, time="03:59:59", memory="6G")
        self._validate_action(action)
        if action in ("pvacseq", "pvacfuse", "pvacsplice"):
            return ResourceUsage(
                threads=min(self.default_resource_usage[action].threads, self.cfg.n_threads),
                time=self.default_resource_usage[action].time,
                memory=self.default_resource_usage[action].memory,
            )
        return self.default_resource_usage[action]


class PvacSeqStepPart(PvacToolsStepPart):
    """
    Specifics for pVACseq:

    - pileups of expression data at somatic variant loci
    - gene- & transcript-based TPM for gene overlapping variant
    """

    #: Step name
    name = "pvacseq"

    #: Actions
    actions = ("pileup", "rename", "combine", "pvacseq")

    #: Resources
    default_resource_usage = {
        "pileup": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "rename": ResourceUsage(
            threads=1,
            time="00:59:59",
            memory="2G",
        ),
        "combine": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "pvacseq": ResourceUsage(
            threads=16,
            time="23:59:59",
            memory="64G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: PVACseqModel = self.config.pvacseq

    def _get_input_files_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {}

        tpl = "output/{mapper}.{library}/out/{mapper}.{library}.bam".format(
            mapper=self.config.pileup.tool_rna_mapping,
            library=self.parent.tumor_rna[wildcards.tumor_dna],
        )
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        input_files["bam"] = ngs_mapping(tpl)

        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.annotated_template)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        input_files["loci"] = annotation(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path

        return input_files

    def _get_input_files_rename(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.annotated_template)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        return {"annotated": annotation(tpl)}

    def _get_input_files_combine(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {
            "annotated": "work/{tpl}/out/{tpl}.renamed.vcf.gz".format(tpl=self.annotated_template)
        }

        if self.config.pileup.enabled and wildcards.tumor_dna in self.parent.tumor_rna:
            tpl = "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.annotated_template)
            input_files["pileup"] = tpl

        if self.config.quantification.enabled and (
            library := self.parent.tumor_rna.get(wildcards.tumor_dna, None)
        ):
            name = f"{self.config.quantification.tool_gene_expression_quantification}.{library}"
            tpl = f"output/{name}/out/{name}.gene.sf"
            quantification = self.parent.sub_workflows["gene_expression_quantification"]
            input_files["gene_tpms"] = quantification(tpl)

            tpl = f"output/{name}/out/{name}.transcript.sf"
            input_files["transcript_tpms"] = quantification(tpl)

        return input_files

    @dictify
    def _get_input_files_pvacseq(self, wildcards: Wildcards):
        if self.cfg.path_container:
            yield "container", self.cfg.path_container
        else:
            yield "container", "work/containers/out/pvactools.sif"

        yield "vcf", "work/{tpl}/out/{tpl}.combined.vcf.gz".format(tpl=self.annotated_template)

        yield "alleles", self._get_hla_files(wildcards)

        if self.config.phasing.enabled:
            yield "phased", "work/{tpl}/out/{tpl}.phased.vcf.gz".format(tpl=self.annotated_template)

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.cfg.run_reference_proteome_similarity and self.cfg.peptide_fasta:
            yield "peptides", self.cfg.peptide_fasta

    def _get_output_files_pileup(self):
        return {"vcf": "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.annotated_template)}

    def _get_output_files_rename(self):
        return {"vcf": "work/{tpl}/out/{tpl}.renamed.vcf.gz".format(tpl=self.annotated_template)}

    def _get_output_files_combine(self):
        return {"vcf": "work/{tpl}/out/{tpl}.combined.vcf.gz".format(tpl=self.annotated_template)}

    def _get_output_files_pvacseq(self):
        return self._get_output_files_run()

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

    def _get_args_rename(self, wildcards: Wildcards) -> dict[str, str]:
        return self._get_sample_names(wildcards)

    def _get_args_combine(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.quantification.model_dump(by_alias=True))

        del args["enabled"]
        del args["path_gene_expression_quantification"]
        args["format"] = args.pop("tool_gene_expression_quantification")

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))

        sample_names = self._get_sample_names(wildcards)
        return {
            "extra_args": extra_args.strip(),
            "tumor_sample": sample_names["tumor_sample"],
            "normal_sample": sample_names["normal_sample"],
        }

    def _get_args_pvacseq(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["peptide_fasta"]
        del args["genes_of_interest_file"]
        if args["net_chop_method"] == NetChopMethod.DISABLED:
            del args["net_chop_method"]
            del args["net_chop_threshold"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

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

        class_i = self._read_hla_values(input["alleles"], MHC_CLASSES[0])
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASSES[1])

        samples = self._get_sample_names(wildcards)

        return {
            "normal_sample": samples["normal_sample"],
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles"],
            "extra_args": extra_args.strip(),
        }


class PvacFuseStepPart(PvacToolsStepPart):
    """Specifics for pVACfuse"""

    #: Step name
    name = "pvacfuse"

    #: Actions
    actions = ("pvacfuse",)

    require_rna: bool = True

    #: Resources
    default_resource_usage = {
        "pvacfuse": ResourceUsage(
            threads=16,
            time="23:59:59",
            memory="64G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: PVACfuseModel = self.config.pvacfuse

    def _get_input_files_pvacfuse(self, wildcards: Wildcards):
        input_files = {}
        if self.cfg.path_container:
            input_files["container"] = self.cfg.path_container
        else:
            input_files["container"] = "work/containers/out/pvactools.sif"

        library = self.parent.tumor_rna.get(wildcards.tumor_dna)
        somatic_gene_fusion_calling = self.parent.sub_workflows["somatic_gene_fusion_calling"]
        tpl = f"{self.cfg.tool_somatic_gene_fusion_calling}.{library}"
        input_files["fusions"] = somatic_gene_fusion_calling(
            "output/" + tpl + "/out/" + tpl + ".fusions.tsv"
        )

        input_files["alleles"] = self._get_hla_files(wildcards)

        if self.cfg.genes_of_interest_file:
            input_files["genes"] = self.cfg.genes_of_interest_file
        if self.cfg.run_reference_proteome_similarity and self.cfg.peptide_fasta:
            input_files["peptides"] = self.cfg.peptide_fasta

        return input_files

    def _get_output_files_pvacfuse(self):
        return self._get_output_files_run()

    def _get_args_pvacfuse(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["path_somatic_gene_fusion_calling"]
        del args["tool_somatic_gene_fusion_calling"]
        del args["peptide_fasta"]
        del args["genes_of_interest_file"]
        if args["net_chop_method"] == NetChopMethod.DISABLED:
            del args["net_chop_method"]
            del args["net_chop_threshold"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

        if args["percentile_threshold"] is None:
            del args["percentile_threshold"]
            del args["percentile_threshold_strategy"]
        else:
            del args["binding_threshold"]

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))

        class_i = self._read_hla_values(input["alleles"], MHC_CLASSES[0])
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASSES[1])

        samples = self._get_sample_names(wildcards)

        return {
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles"],
            "extra_args": extra_args.strip(),
        }


class PvacSpliceStepPart(PvacToolsStepPart):
    """
    Specifics for pVACseq:

    - pileups of expression data at somatic variant loci
    - gene- & transcript-based TPM for gene overlapping variant
    """

    #: Step name
    name = "pvacsplice"

    #: Actions
    actions = ("junction", "pvacsplice")

    require_rna: bool = True

    #: Resources
    default_resource_usage = {
        "junction": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "pvacsplice": ResourceUsage(
            threads=16,
            time="23:59:59",
            memory="64G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: PVACspliceModel = self.config.pvacsplice

    def _get_input_files_junction(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {}

        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.annotated_template)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        input_files["annotated"] = annotation(tpl)

        tpl = "output/{mapper}.{library}/out/{mapper}.{library}.bam".format(
            mapper=self.config.pileup.tool_rna_mapping,
            library=self.parent.tumor_rna[wildcards.tumor_dna],
        )
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        input_files["bam"] = ngs_mapping(tpl)

        tpl = "output/{mapper}.{library}/strandedness/{mapper}.{library}.decision.json".format(
            mapper=self.config.pileup.tool_rna_mapping,
            library=self.parent.tumor_rna[wildcards.tumor_dna],
        )
        input_files["strandedness"] = ngs_mapping(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path
        input_files["features"] = self.w_config.static_data_config.features.path

        return input_files

    @dictify
    def _get_input_files_pvacsplice(self, wildcards: Wildcards):
        if self.cfg.path_container:
            yield "container", self.cfg.path_container
        else:
            yield "container", "work/containers/out/pvactools.sif"

        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

        yield (
            "annotated",
            "work/{tpl}/out/{tpl}.renamed.vcf.gz".format(tpl=self.annotated_template),
        )

        yield (
            "junctions",
            "work/{tpl}/out/{tpl}.junctions.tsv".format(tpl=self.annotated_template),
        )

        yield "alleles", self._get_hla_files(wildcards)

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.cfg.run_reference_proteome_similarity and self.cfg.peptide_fasta:
            yield "peptides", self.cfg.peptide_fasta

    def _get_output_files_junction(self):
        return {
            "junctions": "work/{tpl}/out/{tpl}.junctions.tsv".format(tpl=self.annotated_template)
        }

    def _get_output_files_pvacsplice(self):
        return self._get_output_files_run()

    def _get_args_junction(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        decision = "no file"
        with open(input["strandedness"], "rt") as f:
            decision = json.load(f).get("decision", "not found")
        match decision:
            case "1":
                decision = "FR"
            case "2":
                decision = "RF"
            case _:
                rna_sample = self.parent.tumor_rna.get(wildcards.tumor_dna, wildcards.tumor_dna)
                raise UnsupportedProtocolStrand(
                    f"Illegal strandedness {decision} for sample {rna_sample}"
                )
        return {"strandedness": decision}

    def _get_args_pvacsplice(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["peptide_fasta"]
        del args["genes_of_interest_file"]
        if args["net_chop_method"] == NetChopMethod.DISABLED:
            del args["net_chop_method"]
            del args["net_chop_threshold"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

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

        class_i = self._read_hla_values(input["alleles"], MHC_CLASSES[0])
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASSES[1])

        samples = self._get_sample_names(wildcards)

        return {
            "normal_sample": samples["normal_sample"],
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles"],
            "extra_args": extra_args.strip(),
        }


class PhasingStepPart(BaseStepPart):
    """
    Phase sometic with germline variants using obsolete GATK, for pVACtools only.

    TODO: A better solution sould be developed using current GATK tools
    """

    name = "phasing"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, time="23:59:59", memory="32G")}

    def __init__(self, parent):
        super().__init__(parent)
        self.annotated_template = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            self.annotated_template += ".filtered"
        self.annotated_template += ".{tumor_dna}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path

        combined = self.parent.sub_workflows["combine_variants"]
        tpl = f"{self.config.phasing.tool_ngs_mapping}.combined.{wildcards.tumor_dna}"
        yield "vcf", combined(os.path.join("output", tpl, "out", tpl + ".vcf.gz"))

        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = f"{self.config.phasing.tool_ngs_mapping}.{wildcards.tumor_dna}"
        yield "bam", ngs_mapping(os.path.join("output", tpl, "out", tpl + ".bam"))

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {
                    "vcf": "work/{tpl}/out/{tpl}.phased.vcf.gz".format(tpl=self.annotated_template)
                }
            case _:
                raise MissingConfiguration(f"Unknown action {action} during phasing")

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}

    def get_log_file(self, action):
        """Return mapping of log files."""
        self._validate_action(action)
        tpl = "work/{tpl}/log/{action}".format(tpl=self.annotated_template, action=self.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = tpl + ext
            log_files[key + "_md5"] = log_files[key] + ".md5"
        return log_files

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return self.default_resource_usage[action]


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
        cfg = config["step_config"]["somatic_neoepitope_prediction"]
        if cfg.get("pileup", {}).get("enabled", False) or cfg.get("phasing", {}).get(
            "enabled", False
        ):
            from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

            previous_steps.append(NgsMappingWorkflow)
        if cfg.get("quantification", {}).get("enabled", False):
            from snappy_pipeline.workflows.gene_expression_quantification import (
                GeneExpressionQuantificationWorkflow,
            )

            previous_steps.append(GeneExpressionQuantificationWorkflow)
        if cfg.get("phasing", {}).get("enabled", False):
            from snappy_pipeline.workflows.combine_variants import CombineVariantsWorkflow

            previous_steps.append(CombineVariantsWorkflow)

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticNeoepitopePredictionConfigModel,
            previous_steps=previous_steps,
        )

        self.register_sub_step_classes(
            (
                PvacToolsStepPart,
                PvacSeqStepPart,
                PvacFuseStepPart,
                PvacSpliceStepPart,
                PhasingStepPart,
                LinkOutStepPart,
            )
        )

        self.sample_table: pd.DataFrame = sample_sheets(self.sheets)
        assert "extractionType" in self.sample_table.columns, (
            "'extractionType' missing from sample sheet"
        )
        self.tumor_dna = tumor_to_normal_mapping(
            self.sample_table[self.sample_table["extractionType"] == ExtractionType.DNA]
        )
        self.tumor_rna = self._dna_to_rna_mapping(self.sample_table)
        if (self.config.pileup.enabled or self.config.quantification.enabled) or (
            "pvacfuse" in self.config.tools or "pvacsplice" in self.config.tools
        ):
            assert any(map(lambda lib: lib in self.tumor_rna.keys(), self.tumor_dna.keys())), (
                "No tumor sample with somatic variant has expression data"
            )

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
        if self.config.phasing.enabled:
            self.register_sub_workflow(
                "combine_variants",
                self.config.phasing.path_combine_variants,
            )
        if "pvacfuse" in self.config.tools:
            self.register_sub_workflow(
                "somatic_gene_fusion_calling",
                self.config.pvacfuse.path_somatic_gene_fusion_calling,
            )

    def get_result_files(self):
        all_tools = []
        for tool in self.config.tools:
            all_tools += self.sub_steps[tool].get_result_files()
        return all_tools

    def check_config(self):
        hla_typing_config = self.w_config.step_config.get("hla_typing", None)
        for extraction_type in EXTRACTION_TYPES:
            for mhc_class in MHC_CLASSES:
                tool = self.config.tools_hla_typing.get(extraction_type, {}).get(mhc_class, None)
                if tool and tool not in hla_typing_config.tools.get(extraction_type, []):
                    raise MissingConfiguration(f"hla_typing tool {tool} not configured")

    def _dna_to_rna_mapping(self, sample_table: pd.DataFrame) -> dict[str, str]:
        dna = sample_table[sample_table["extractionType"] == ExtractionType.DNA]
        rna = sample_table[sample_table["extractionType"] == ExtractionType.RNA]
        dna_rna_map = dna[["ngs_library", "bio_entity", "bio_sample"]].merge(
            rna[["ngs_library", "bio_entity", "bio_sample"]], on=["bio_entity", "bio_sample"]
        )
        return pd.Series(
            dna_rna_map.ngs_library_y.values, index=dna_rna_map.ngs_library_x.values
        ).to_dict()

# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_neoepitope_prediction`` step

The somatic_neoepitope_prediction step allows for the prediction of neoepitopes from somatic
(small) variant calling results, from splicing events and from gene fusions. The latter two
require that RNA data is present beside tumor/normal pair.

The current implementation is using pVACtools's "pVACseq", "pVACsplice" & "pVACfuse", _via_
its docker image. The IEDB tools within the image are *NOT* free software, and their use
is restricted to academic, non-profit research (unless a paying licence is acquired).

The netchop & netMHCstab tools are *NOT* included in the docker image, and must be run separately.
(the docker image *CAN* will these tools, but by making remote calls to the DTU servers).
We have an experimental implementation which can use a local installation of "NetChop",
but it is largely untested. In particular, it is unclear if it considers protein sequence changes
due to germline variants in the vicinity of somatic ones.

.. note::

    Status: under development

==========
Step Input
==========

Multiple steps are required to benefit from all features offered by pVACtools. Each tools has different
requirements.

-------
pVACseq
-------

- `hla_typing` (required): HLA types (both MHC classes I & II).
- `somatic_variant_annotation` (required): somatic variants.
- `gene_expression_quantification` (optional): expression TPMs computed by the `salmon` tool.
- `ngs_mapping` (optional): somatic mutation experimental evidence as read counts from `star` tool.
- `combine_variants` (optional): allow to take into account protein sequence changes due to germline variants in the vicinity of somatic ones.
- `create_proteome` (optional, not yet connected): create personalised proteome to search neoepitope sequences.

----------
pVACsplice
----------

- `hla_typing` (required): HLA types (both MHC classes I & II).
- `somatic_variant_annotation` (required): somatic variants.
- `ngs_mapping` (required): splicing & regulatory events from `star` & `strandedness` tools.
- `create_proteome` (optional, not yet connected): create personalised proteome to search neoepitope sequences.

--------
pVACfuse
--------

- `hla_typing` (required): HLA types (both MHC classes I & II).
- `somatic_gene_fusion_calling` (required): gene fusions found in RNA data.
- `create_proteome` (optional, not yet connected): create personalised proteome to search neoepitope sequences.

Note that germline variants cannot be obtained using the `variant_calling` step, because of biomedsheet incompatibilites.
They need be be called using `germline_variant_calling`. Note also that it is advisable to filter
both somatic & germline variants before they are combined.

===========
Step Output
===========

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_neoepitope_prediction.rst

"""

import json
import os

import pandas as pd

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.iocontainers import Wildcards, InputFiles

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.models.common import ExtractionType
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.variant_annotation import VariantAnnotationWorkflow
from snappy_pipeline.workflows.variant_calling import VariantCallingWorkflow
from snappy_pipeline.workflows.variant_filtration import VariantFiltrationWorkflow
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.gene_expression_quantification import (
    GeneExpressionQuantificationWorkflow,
)
from snappy_pipeline.workflows.combine_variants import CombineVariantsWorkflow
from snappy_pipeline.workflows.somatic_gene_fusion_calling import SomaticGeneFusionCallingWorkflow
from .model import (
    SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel,
    SupportedPredictionTool,
)
from .model import PVACseq as PVACseqModel
from .model import PVACfuse as PVACfuseModel
from .model import PVACsplice as PVACspliceModel
from .model import NetChop as NetChopModel
from .model import MHC_CLASS, MHC_CLASS_I, MHC_CLASS_II


__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

#: Extensions of files to create as main payload
PREPARE_EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticNeoepitopePredictionConfigModel.default_config_yaml_string()


class UnsupportedProtocolStrand(Exception):
    pass


class PvacToolsStepPart(BaseStepPart):
    """Generic stuff for pVACtools modules"""

    #: Step name
    name = "pvactools"

    actions = ("install", "normalize", "normalize_full")
    require_rna: bool = False

    default_resource_usage = {
        "install": ResourceUsage(threads=1, runtime="03:59:59", mem="64G"),
        "normalize": ResourceUsage(threads=1, runtime="01:00:00", mem="4G"),
        "normalize_full": ResourceUsage(threads=1, runtime="01:00:00", mem="4G"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.prepare_tpl = "{tumor_dna}"
        self.output_tpl = f"{{tumor_dna}}.{self.name}"

        self.hla_tools = {}
        for extraction_type in ExtractionType:
            extraction_type = extraction_type.lower()
            for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
                if tool := self.config.tool_hla_typing.get(extraction_type, {}).get(
                    mhc_class.name, None
                ):
                    if extraction_type not in self.hla_tools:
                        self.hla_tools[extraction_type] = {}
                    if tool:
                        try:
                            hla_config = self.parent.get_task_config("hla_typing")
                            mapper = getattr(hla_config, "mapper", None)
                        except Exception:
                            mapper = None
                        if mapper:
                            self.hla_tools[extraction_type][mhc_class.name] = f"{mapper}.{tool}"
                        else:
                            self.hla_tools[extraction_type][mhc_class.name] = tool

        if self.config.proteome.enabled:
            if self.config.proteome.add_unmutated:
                self.proteome_file = "work/pvactools/out/proteome.fa.gz"
            else:
                self.proteome_file = self.config.proteome.external_proteome
        else:
            self.proteome_file = None

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _get_input_files_normalize(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.upstream("somatic_variant_annotation")
        return {"annotated": annotation(tpl)}

    def _get_input_files_normalize_full(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "output/{tpl}/out/{tpl}.full.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.upstream("somatic_variant_annotation")
        return {"annotated": annotation(tpl)}

    def get_output_files(self, action):
        if action == "install":
            return {"container": "work/containers/out/pvactools.sif"}
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def _get_output_files_run(self):
        return {
            "filtered": "work/"
            + self.output_tpl
            + "/out/{mhc_class_d,MHC_Class_II?|combined}/{tumor_dna}.{mhc_class_fn,MHC_II?|Combined}.filtered.tsv",
            "done": "work/"
            + self.output_tpl
            + "/out/{mhc_class_d,MHC_Class_II?|combined}.{tumor_dna}.{mhc_class_fn,MHC_II?|Combined}.done",
        }

    def _get_output_files_normalize(self):
        tpl = "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl)
        return {"vcf": tpl}

    def _get_output_files_normalize_full(self):
        tpl = "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl)
        return {"vcf": tpl}

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_normalize(self, wildcards: Wildcards) -> dict[str, str]:
        return self._get_sample_names(wildcards)

    def _get_args_normalize_full(self, wildcards: Wildcards) -> dict[str, str]:
        return self._get_args_normalize(wildcards)

    def get_log_file(self, action):
        """Return mapping of log files."""
        if action == "install":
            return f"work/containers/log/{self.name}.log"
        if action == self.name:
            tpl = "work/{tpl}/log/{{mhc_class_d,MHC_Class_II?|combined}}.{{tumor_dna}}.{{mhc_class_fn,MHC_II?|Combined}}.filtered".format(
                tpl=self.output_tpl
            )
        else:
            self._validate_action(action)
            tpl = "work/{tpl}/log/{action}".format(tpl=self.prepare_tpl, action=action)
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

    def _get_sample_names(self, wildcards: Wildcards) -> dict[str, str]:
        args = {}
        args["tumor_sample"] = wildcards.tumor_dna
        args["tumor_library"] = wildcards.tumor_dna
        if normal_dna := self.parent.tumor_dna.get(wildcards.tumor_dna, None):
            args["normal_sample"] = normal_dna
            args["normal_library"] = normal_dna
        return args

    @listify
    def _get_hla_files(self, wildcards: Wildcards):
        hla_typing = self.parent.upstream("hla_typing")
        tumor_dna = wildcards.tumor_dna
        normal_dna = self.parent.tumor_dna.get(tumor_dna, None)
        tumor_rna = self.parent.tumor_rna.get(tumor_dna, None)

        input_files = []
        for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
            if tool := self.hla_tools.get("dna", {}).get(mhc_class.name, None):
                tpl = "output/{tool}.{library_name}/out/{tool}.{library_name}.json"
                input_files.append(tpl.format(tool=tool, library_name=tumor_dna))
                if normal_dna:
                    input_files.append(tpl.format(tool=tool, library_name=normal_dna))
        if tumor_rna:
            for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
                if tool := self.hla_tools.get("rna", {}).get(mhc_class.name, None):
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
    def _read_hla_values(cls, hla_typing_files: list[str], mhc_class: MHC_CLASS) -> list[str]:
        hla_types = []
        for fn in hla_typing_files:
            with open(fn, "rt") as f:
                calls: dict[str, Any] = json.load(f)
            for locus, alleles in calls.items():
                if locus in mhc_class.genes:
                    for allele in alleles:
                        m = mhc_class.pattern.match(allele)
                        if m:
                            hla_types.append(mhc_class.prefix + m.group("valid"))

        return list(set(hla_types))

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        if action in ("pvacseq", "pvacfuse", "pvacsplice"):
            return ResourceUsage(
                threads=min(self.default_resource_usage[action].threads, self.cfg.n_threads),
                runtime=self.default_resource_usage[action].runtime,
                mem=self.default_resource_usage[action].mem,
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
    actions = ("pileup", "combine", "pvacseq")

    #: Resources
    default_resource_usage = {
        "pileup": ResourceUsage(
            threads=1,
            runtime="03:59:59",
            mem="6G",
        ),
        "combine": ResourceUsage(
            threads=1,
            runtime="03:59:59",
            mem="6G",
        ),
        "pvacseq": ResourceUsage(
            threads=16,
            runtime="23:59:59",
            mem="64G",
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
        ngs_mapping = self.parent.upstream("ngs_mapping")
        input_files["bam"] = ngs_mapping(tpl)

        if self.cfg.use_all_transcripts:
            tpl = "output/{tpl}/out/{tpl}.full.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.upstream("somatic_variant_annotation")
        input_files["loci"] = annotation(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path

        return input_files

    def _get_input_files_combine(self, wildcards: Wildcards) -> dict[str, str]:
        if self.cfg.use_all_transcripts:
            tpl = "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            tpl = "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl)
        input_files = {"annotated": tpl}

        if self.config.pileup.enabled and wildcards.tumor_dna in self.parent.tumor_rna:
            tpl = "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.prepare_tpl)
            input_files["pileup"] = tpl

        if self.config.quantification.enabled and (
            library := self.parent.tumor_rna.get(wildcards.tumor_dna, None)
        ):
            name = f"{self.config.quantification.tool_gene_expression_quantification}.{library}"
            tpl = f"output/{name}/out/{name}.gene.sf"
            quantification = self.parent.upstream("gene_expression_quantification")
            input_files["gene_tpms"] = quantification(tpl)

            tpl = f"output/{name}/out/{name}.transcript.sf"
            input_files["transcript_tpms"] = quantification(tpl)

            if self.config.quantification.duplicate_transcripts_table:
                input_files["duplicate_transcripts_table"] = (
                    self.config.quantification.duplicate_transcripts_table
                )

        return input_files

    @dictify
    def _get_input_files_pvacseq(self, wildcards: Wildcards):
        if self.cfg.path_container:
            yield "container", self.cfg.path_container
        else:
            yield "container", "work/containers/out/pvactools.sif"

        if wildcards.tumor_dna in self.parent.tumor_rna and (
            self.config.pileup.enabled or self.config.quantification.enabled
        ):
            yield "vcf", "work/{tpl}/out/{tpl}.combined.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            if self.cfg.use_all_transcripts:
                yield (
                    "vcf",
                    "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl),
                )
            else:
                yield "vcf", "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl)

        yield "alleles", self._get_hla_files(wildcards)

        if self.config.phasing.enabled:
            yield "phased", "work/{tpl}/out/{tpl}.phased.vcf.gz".format(tpl=self.prepare_tpl)

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.proteome_file:
            yield "peptides", self.proteome_file

    def _get_output_files_pileup(self):
        return {"vcf": "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.prepare_tpl)}

    def _get_output_files_combine(self):
        return {"vcf": "work/{tpl}/out/{tpl}.combined.vcf.gz".format(tpl=self.prepare_tpl)}

    def _get_output_files_pvacseq(self):
        return self._get_output_files_run()

    def _get_args_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.pileup.model_dump(by_alias=True))

        del args["enabled"]
        del args["path_ngs_mapping"]
        del args["tool_rna_mapping"]
        if args["baq"] is None:
            del args["baq"]

        args = self._extra_args_lists(args)
        extra_args = " ".join(sorted(list(self._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(self._group_extra_args(args))))

        return {
            "tumor_sample": self._get_sample_names(wildcards)["tumor_sample"],
            "extra_args": extra_args.strip(),
        }

    def _get_args_combine(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.quantification.model_dump(by_alias=True))

        del args["enabled"]
        del args["path_gene_expression_quantification"]
        del args["duplicate_transcripts_table"]
        args["format"] = args.pop("tool_gene_expression_quantification")

        args = self._extra_args_lists(args)
        extra_args = " ".join(sorted(list(self._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(self._group_extra_args(args))))

        sample_names = self._get_sample_names(wildcards)
        return {
            "extra_args": extra_args.strip(),
            "tumor_sample": sample_names["tumor_sample"],
            "normal_sample": sample_names["normal_sample"],
        }

    def _get_args_pvacseq(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["use_all_transcripts"]
        del args["genes_of_interest_file"]
        del args["net_chop"]
        del args["netmhc_stab"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

        if args["maximum_transcript_support_level"] is None:
            del args["maximum_transcript_support_level"]

        samples = self._get_sample_names(wildcards)

        if "normal_sample" not in samples:
            del args["normal_vaf"]
            del args["normal_cov"]

        args = self._extra_args_lists(args)
        extra_args = " ".join(sorted(list(self._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(self._group_extra_args(args))))

        class_i = self._read_hla_values(input["alleles"], MHC_CLASS_I)
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASS_II)
        return {
            "normal_sample": samples.get("normal_sample", None),
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles", "filtered"],
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
            runtime="23:59:59",
            mem="64G",
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
        somatic_gene_fusion_calling = self.parent.upstream("somatic_gene_fusion_calling")
        tpl = f"{self.cfg.tool_somatic_gene_fusion_calling}.{library}"
        input_files["fusions"] = somatic_gene_fusion_calling(
            "output/" + tpl + "/out/" + tpl + ".fusions.tsv"
        )

        input_files["alleles"] = self._get_hla_files(wildcards)

        if self.cfg.genes_of_interest_file:
            input_files["genes"] = self.cfg.genes_of_interest_file
        if self.proteome_file:
            input_files["peptides"] = self.proteome_file

        return input_files

    def _get_output_files_pvacfuse(self):
        return self._get_output_files_run()

    def _get_args_pvacfuse(self, wildcards: Wildcards, input: InputFiles) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        del args["path_container"]
        del args["path_somatic_gene_fusion_calling"]
        del args["tool_somatic_gene_fusion_calling"]
        del args["net_chop"]
        del args["netmhc_stab"]
        del args["genes_of_interest_file"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

        args = self._extra_args_lists(args)
        extra_args = " ".join(sorted(list(self._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(self._group_extra_args(args))))

        class_i = self._read_hla_values(input["alleles"], MHC_CLASS_I)
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASS_II)

        samples = self._get_sample_names(wildcards)

        return {
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles", "filtered"],
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
    actions = ("junction", "pvacsplice", "workaround")

    require_rna: bool = True

    #: Resources
    default_resource_usage = {
        "junction": ResourceUsage(
            threads=1,
            runtime="03:59:59",
            mem="16G",
        ),
        "pvacsplice": ResourceUsage(
            threads=16,
            runtime="23:59:59",
            mem="64G",
        ),
        "workaround": ResourceUsage(
            threads=1,
            runtime="01:00:00",
            mem="4G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: PVACspliceModel = self.config.pvacsplice

    def _get_input_files_junction(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {}

        if self.cfg.use_all_transcripts:
            tpl = "output/{tpl}/out/{tpl}.full.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.upstream("somatic_variant_annotation")
        input_files["annotated"] = annotation(tpl)

        tpl = "output/{mapper}.{library}/out/{mapper}.{library}.bam".format(
            mapper=self.config.pileup.tool_rna_mapping,
            library=self.parent.tumor_rna[wildcards.tumor_dna],
        )
        ngs_mapping = self.parent.upstream("ngs_mapping")
        input_files["bam"] = ngs_mapping(tpl)

        tpl = "output/{mapper}.{library}/strandedness/{mapper}.{library}.decision.json".format(
            mapper=self.config.pileup.tool_rna_mapping,
            library=self.parent.tumor_rna[wildcards.tumor_dna],
        )
        input_files["strandedness"] = ngs_mapping(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path
        # input_files["features"] = self.w_config.static_data_config.features.path
        input_files["features"] = "work/pvacsplice_workaround/out/features.gtf"

        return input_files

    @dictify
    def _get_input_files_pvacsplice(self, wildcards: Wildcards):
        if self.cfg.path_container:
            yield "container", self.cfg.path_container
        else:
            yield "container", "work/containers/out/pvactools.sif"

        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

        if self.cfg.use_all_transcripts:
            yield (
                "annotated",
                "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl),
            )
        else:
            yield (
                "annotated",
                "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl),
            )

        yield (
            "junctions",
            "work/{tpl}/out/{tpl}.junctions.tsv".format(tpl=self.prepare_tpl),
        )

        yield "alleles", self._get_hla_files(wildcards)

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.proteome_file:
            yield "peptides", self.proteome_file

    def _get_output_files_junction(self):
        return {"junctions": "work/{tpl}/out/{tpl}.junctions.tsv".format(tpl=self.prepare_tpl)}

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
        del args["use_all_transcripts"]
        del args["net_chop"]
        del args["netmhc_stab"]
        del args["genes_of_interest_file"]
        n_threads = args.pop("n_threads")

        algorithms = args.pop("algorithms")
        if isinstance(algorithms, list):
            algorithms = " ".join(algorithms)

        if args["maximum_transcript_support_level"] is None:
            del args["maximum_transcript_support_level"]

        args = PvacSeqStepPart._extra_args_lists(args)
        extra_args = " ".join(sorted(list(PvacSeqStepPart._extra_args_flags(args))))

        extra_args += " " + " ".join(sorted(list(PvacSeqStepPart._group_extra_args(args))))

        class_i = self._read_hla_values(input["alleles"], MHC_CLASS_I)
        class_ii = self._read_hla_values(input["alleles"], MHC_CLASS_II)

        samples = self._get_sample_names(wildcards)

        return {
            "normal_sample": samples["normal_sample"],
            "tumor_sample": samples["tumor_sample"],
            "class_i": sorted(class_i),
            "class_ii": sorted(class_ii),
            "algorithms": algorithms,
            "n_threads": n_threads,
            "exclude_bind": ["container", "alleles", "filtered"],
            "extra_args": extra_args.strip(),
        }


class PhasingStepPart(BaseStepPart):
    """
    Phase somatic with germline variants using obsolete GATK, for pVACtools only.

    TODO: A better solution sould be developed using current GATK tools
    """

    name = "phasing"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, runtime="23:59:59", mem="32G")}

    def __init__(self, parent):
        super().__init__(parent)
        self.prepare_tpl = "{tumor_dna}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "reference", self.w_config.static_data_config.reference.path

        combined = self.parent.upstream("combine_variants")
        tpl = f"{self.config.phasing.tool_ngs_mapping}.combined.{wildcards.tumor_dna}"
        yield "vcf", combined(os.path.join("output", tpl, "out", tpl + ".vcf.gz"))

        ngs_mapping = self.parent.upstream("ngs_mapping")
        tpl = f"{self.config.phasing.tool_ngs_mapping}.{wildcards.tumor_dna}"
        yield "bam", ngs_mapping(os.path.join("output", tpl, "out", tpl + ".bam"))

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {"vcf": "work/{tpl}/out/{tpl}.phased.vcf.gz".format(tpl=self.prepare_tpl)}
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
        tpl = "work/{tpl}/log/{action}".format(tpl=self.prepare_tpl, action=self.name)
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


class NetChopStepPart(BaseStepPart):
    """run netchop after pVACseq, pVACsplice or pVACfuse"""

    name = "netchop"
    actions = ("pvacseq", "pvacfuse", "pvacsplice")

    #: Resources
    resource_usage = ResourceUsage(threads=1, runtime="23:59:59", mem="32G")
    default_resource_usage = {
        "pvacseq": resource_usage,
        "pvacfuse": resource_usage,
        "pvacsplice": resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.prepare_tpl = "{tumor_dna}"
        self.input_tpl = "{tool}.{tumor_dna}"
        self.output_tpl = "{tool}.{tumor_dna}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _get_input_files_pvacseq(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacseq", (
            "Internal error: tool {wildcards.tool} should be 'pvacseq'"
        )
        return self._get_input_files_run(wildcards)

    def _get_input_files_pvacfuse(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacfuse", (
            "Internal error: tool {wildcards.tool} should be 'pvacfuse'"
        )
        return self._get_input_files_run(wildcards)

    def _get_input_files_pvacsplice(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacsplice", (
            "Internal error: tool {wildcards.tool} should be 'pvacsplice'"
        )
        return self._get_input_files_run(wildcards)

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        yield "netchop", self.config.get(wildcards.tool).net_chop.path_netchop

        if self.config.phasing.enabled:
            combined = self.parent.upstream("combine_variants")
            tpl = f"{self.config.phasing.tool_ngs_mapping}.combined.{wildcards.tumor_dna}"
            yield "vcf", combined(os.path.join("output", tpl, "out", tpl + ".vcf.gz"))
        else:
            if self.config.get(wildcards.tool).use_all_transcripts:
                yield (
                    "vcf",
                    "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl),
                )
            else:
                yield "vcf", "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl)

        tool = wildcards.tool
        if self.config.get(tool).class_i_epitope_length:
            tool_dirname = "MHC_Class_I"
            tool_filename = "MHC_I"
            if self.config.get(tool).class_ii_epitope_length:
                tool_dirname = "combined"
                tool_filename = "Combined"
        else:
            tool_dirname = "MHC_Class_II"
            tool_filename = "MHC_II"

        yield (
            "epitopes",
            os.path.join(
                "work",
                self.input_tpl.format(**wildcards),
                "out",
                tool_dirname,
                f"{wildcards.tumor_dna}.{tool_filename}.filtered.tsv",
            ),
        )

    def get_output_files(self, action: str) -> dict[str, Any]:
        match action:
            case "pvacseq" | "pvacfuse" | "pvacsplice":
                return {
                    "epitopes": f"work/{self.output_tpl}/out/{{mhc_class_d}}/{{sample}}.{{mhc_class_fn}}.netchop.tsv"
                }
            case _:
                raise MissingConfiguration(f"Unknown action {action} during phasing")

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_pvacseq(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacseq", (
            "Internal error: tool {wildcards.tool} should be 'pvacseq'"
        )
        return self._get_args_run(wildcards)

    def _get_args_pvacfuse(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacfuse", (
            "Internal error: tool {wildcards.tool} should be 'pvacfuse'"
        )
        return self._get_args_run(wildcards)

    def _get_args_pvacsplice(self, wildcards: Wildcards):
        assert wildcards.tool == "pvacsplice", (
            "Internal error: tool {wildcards.tool} should be 'pvacsplice'"
        )
        return self._get_args_run(wildcards)

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg: NetChopModel = self.config.get(wildcards.tool).net_chop
        return {"tool": wildcards.tool, "method": cfg.method, "threshold": cfg.threshold}

    def get_log_file(self, action):
        """Return mapping of log files."""
        self._validate_action(action)
        tpl = f"work/{self.output_tpl}/log/{{mhc_class_d}}.{{sample}}.{{mhc_class_fn}}.netchop"
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


class ProteomeStepPart(BaseStepPart):
    """Create personalised proteome"""

    name = "proteome"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, runtime="03:59:59", mem="24G")}

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        cfg = self.config.proteome
        if cfg.external_proteome:
            yield "path_proteome", cfg.external_proteome
        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

    def get_output_files(self, action: str) -> dict[str, str]:
        match action:
            case "run":
                return {"proteome": "work/pvactools/out/proteome.fa.gz"}
            case _:
                raise MissingConfiguration(
                    f"Unknown action {action} during personal proteome building"
                )

    def get_args(self, action: str) -> dict[str, Any]:
        match action:
            case "run":
                return {"add_unmutated": self.config.proteome.add_unmutated}
            case _:
                raise MissingConfiguration(
                    f"Unknown action {action} during personal proteome building"
                )

    def get_log_file(self, action: str) -> dict[str, str]:
        self._validate_action(action)
        tpl = "work/pvactools/log/proteome"
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

    produces = [DataSignature(DataType.TABULAR, frozenset({"neoepitope"}))]
    config_model_class = SomaticNeoepitopePredictionConfigModel

    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir, **kwargs):
        previous_steps = [
            VariantAnnotationWorkflow,
            HlaTypingWorkflow,
            VariantCallingWorkflow,
            VariantFiltrationWorkflow,
            NgsMappingWorkflow,
            GeneExpressionQuantificationWorkflow,
            CombineVariantsWorkflow,
            SomaticGeneFusionCallingWorkflow,
        ]

        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            previous_steps=previous_steps,
            **kwargs,
        )

        match self.config.tool:
            case SupportedPredictionTool.PVACSEQ:
                selected_sub_step = PvacSeqStepPart
            case SupportedPredictionTool.PVACFUSE:
                selected_sub_step = PvacFuseStepPart
            case SupportedPredictionTool.PVACSPLICE:
                selected_sub_step = PvacSpliceStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes(
            (
                PvacToolsStepPart,
                selected_sub_step,
                PhasingStepPart,
                NetChopStepPart,
                ProteomeStepPart,
                LinkOutStepPart,
            )
        )

        df = self.build_library_dataframe()
        assert "extraction_type" in df.columns, "'extraction_type' missing from library dataframe"

        tumor_df = df[df["role"] == "tumor"]
        self.tumor_dna = dict(zip(tumor_df["library_name"], tumor_df.get("matched_normal_lib", "")))
        self.tumor_dna = {k: v for k, v in self.tumor_dna.items() if v}

        self.tumor_rna = self._dna_to_rna_mapping(df)

        if (self.config.pileup.enabled or self.config.quantification.enabled) or (
            self.config.tool
            in (SupportedPredictionTool.PVACFUSE, SupportedPredictionTool.PVACSPLICE)
        ):
            assert any(lib in self.tumor_rna for lib in self.tumor_dna), (
                "No tumor sample with somatic variant has expression data"
            )

    @listify
    def get_result_files(self):
        log_exts = ("log", "conda_list.txt", "conda_info.txt")
        hash_exts = ("", ".md5")

        tumor_samples = [lib for lib, norm in self.tumor_dna.items() if norm]

        tool_name = self.config.tool
        tool = self.sub_steps[tool_name]
        tool_cfg = self.config.get(tool_name)

        if tool_cfg.class_i_epitope_length:
            mhc_class_d = "MHC_Class_I"
            mhc_class_fn = "MHC_I"
            if tool_cfg.class_ii_epitope_length:
                mhc_class_d = "combined"
                mhc_class_fn = "Combined"
        else:
            mhc_class_d = "MHC_Class_II"
            mhc_class_fn = "MHC_II"

        if tool_cfg.net_chop.enabled:
            ext = "netchop"
        else:
            ext = "filtered"

        for tumor_dna in tumor_samples:
            if tool.require_rna and self.tumor_rna.get(tumor_dna, None) is None:
                continue

            d = f"output/{tool_name}.{tumor_dna}"
            fn = f"out/{mhc_class_d}/{tumor_dna}.{mhc_class_fn}.{ext}.tsv"
            yield f"{d}/{fn}"

            fn = f"log/{mhc_class_d}.{tumor_dna}.{mhc_class_fn}.{ext}.{{log_ext}}{{hash_ext}}"
            for log_ext in log_exts:
                for hash_suffix in hash_exts:
                    yield f"{d}/{fn.format(log_ext=log_ext, hash_ext=hash_suffix)}"

    def check_config(self):
        for extraction_type in (ExtractionType.DNA, ExtractionType.RNA):
            for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
                tool = self.config.tool_hla_typing.get(extraction_type, {}).get(
                    mhc_class.name, None
                )
                if tool:
                    self.ensure_w_config(
                        ("static_data_config", "reference", "path"),
                        "Path to reference FASTA not configured but required for neoepitope prediction",
                    )

    def _dna_to_rna_mapping(self, sample_table: pd.DataFrame) -> dict[str, str]:
        dna = sample_table[sample_table["extraction_type"] == ExtractionType.DNA]
        rna = sample_table[sample_table["extraction_type"] == ExtractionType.RNA]
        dna_rna_map = dna[["library_name", "donor_name"]].merge(
            rna[["library_name", "donor_name"]], on=["donor_name"]
        )
        return pd.Series(
            dna_rna_map.library_name_y.values, index=dna_rna_map.library_name_x.values
        ).to_dict()

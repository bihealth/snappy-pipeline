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
from snappy_pipeline.workflows.abstract.exceptions import InvalidConfigurationException
from snappy_pipeline.workflows.common.samplesheet import sample_sheets, tumor_to_normal_mapping
from snappy_pipeline.workflows.somatic_variant_annotation import SomaticVariantAnnotationWorkflow
from snappy_pipeline.workflows.hla_typing import HlaTypingWorkflow
from .model import SomaticNeoepitopePrediction as SomaticNeoepitopePredictionConfigModel
from .model import PVACseq as PVACseqModel
from .model import PVACfuse as PVACfuseModel
from .model import PVACsplice as PVACspliceModel
from .model import NetChop as NetChopModel
from .model import Proteome as ProteomeModel
from .model import GermlineVariantStep
from .model import MHC_CLASS, MHC_CLASS_I, MHC_CLASS_II, ClassIAlgorithm, ClassIIAlgorithm


__author__ = "Eric Blanc"
__email__ = "eric.blanc@bih-charite.de"

#: Extensions of files to create as main payload
PREPARE_EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")
#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticNeoepitopePredictionConfigModel.default_config_yaml_string()


class UnsupportedProtocolStrand(Exception):
    pass


class GenericNeoepitopeStepPart(BaseStepPart):
    # Expressed genes (not genes fragments, not pseudogenes) taken from
    # Marsh et al. (2026) "Nomenclature for Factors of the HLA System, 2026". HLA 107(3):e70595
    # https://doi.org/10.1111/tan.70595

    def __init__(self, parent):
        super().__init__(parent)

    def get_sample_names(self, wildcards: Wildcards) -> dict[str, str]:
        args = {}
        args["tumor_sample"] = wildcards.tumor_dna
        args["tumor_library"] = wildcards.tumor_dna
        if normal_dna := self.parent.tumor_dna.get(wildcards.tumor_dna, None):
            args["normal_sample"] = normal_dna
            args["normal_library"] = normal_dna
        return args


class PvacToolsStepPart(GenericNeoepitopeStepPart):
    """Generic stuff for pVACtools modules"""

    #: Step name
    name = "pvactools"

    actions = ("install", "normalize", "normalize_full")
    require_rna: bool = False

    default_resource_usage = {
        "install": ResourceUsage(threads=1, time="03:59:59", memory="64G"),
        "normalize": ResourceUsage(threads=1, time="01:00:00", memory="4G"),
        "normalize_full": ResourceUsage(threads=1, time="01:00:00", memory="4G"),
    }

    SUBDIRECTORIES = (
        ("MHC_Class_I", "MHC_I"),
        ("MHC_Class_II", "MHC_II"),
        ("combined", "Combined"),
    )
    FILE_EXTENSIONS = (
        ("all", "all_epitopes.tsv"),
        ("filtered", "filtered.tsv"),
        ("aggregated", "all_epitopes.aggregated.tsv"),
        ("json", "all_epitopes.aggregated.metrics.json"),
    )

    def __init__(self, parent):
        super().__init__(parent)
        prefix = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            postfix = "filtered.{tumor_dna}"
        else:
            postfix = "{tumor_dna}"
        self.prepare_tpl = f"{prefix}.{postfix}"
        self.output_tpl = f"{prefix}.{self.name}.{postfix}"

        if self.config.proteome.enabled:
            if self.config.proteome.path_germline_variants:
                self.proteome_file = (
                    f"work/{self.prepare_tpl}/out/{self.prepare_tpl}.proteome.fa.gz"
                )
            elif self.config.proteome.add_unmutated:
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
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        return {"annotated": annotation(tpl)}

    def _get_input_files_normalize_full(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "output/{tpl}/out/{tpl}.full.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        return {"annotated": annotation(tpl)}

    def get_output_files(self, action):
        if action == "install":
            return {"container": "work/containers/out/pvactools.sif"}
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    def _get_output_files_run(self):
        done = ".done"
        outputs = {"done": os.path.join("work", self.output_tpl, "out", done)}
        for k, ext in self.FILE_EXTENSIONS:
            if self.name != "pvacseq" and k == "json":
                continue
            for mhc_class_d, mhc_class_fn in self.SUBDIRECTORIES:
                outputs[f"{k}.{mhc_class_fn}"] = os.path.join(
                    "work",
                    self.output_tpl,
                    "out",
                    mhc_class_d,
                    "{tumor_dna}." + mhc_class_fn + "." + ext,
                )
        return outputs

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
        return self.get_sample_names(wildcards)

    def _get_args_normalize_full(self, wildcards: Wildcards) -> dict[str, str]:
        return self._get_args_normalize(wildcards)

    def get_log_file(self, action):
        """Return mapping of log files."""
        if action == "install":
            return f"work/containers/log/{self.name}.log"
        if action == self.name:
            tpl = f"work/{self.output_tpl}/log/{self.name}.{{tumor_dna}}"
        else:
            self._validate_action(action)
            tpl = "work/{tpl}/log/{action}.{{tumor_dna}}".format(
                tpl=self.prepare_tpl, action=action
            )
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
    actions = ("pileup", "add_expression", "pvacseq")

    #: Resources
    default_resource_usage = {
        "pileup": ResourceUsage(
            threads=1,
            time="03:59:59",
            memory="6G",
        ),
        "add_expression": ResourceUsage(
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

        if self.cfg.use_all_transcripts:
            tpl = "output/{tpl}/out/{tpl}.full.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            tpl = "output/{tpl}/out/{tpl}.vcf.gz".format(tpl=self.prepare_tpl)
        annotation = self.parent.sub_workflows["somatic_variant_annotation"]
        input_files["loci"] = annotation(tpl)

        input_files["reference"] = self.w_config.static_data_config.reference.path

        return input_files

    def _get_input_files_add_expression(self, wildcards: Wildcards) -> dict[str, str]:
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
            quantification = self.parent.sub_workflows["gene_expression_quantification"]
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
            yield "vcf", "work/{tpl}/out/{tpl}.with_expression.vcf.gz".format(tpl=self.prepare_tpl)
        else:
            if self.cfg.use_all_transcripts:
                yield (
                    "vcf",
                    "work/{tpl}/out/{tpl}.normalized.full.vcf.gz".format(tpl=self.prepare_tpl),
                )
            else:
                yield "vcf", "work/{tpl}/out/{tpl}.normalized.vcf.gz".format(tpl=self.prepare_tpl)

        yield "alleles", f"work/{self.prepare_tpl}/out/{self.prepare_tpl}.hla_types.txt"

        if self.config.phasing.enabled:
            yield "phased", "work/{tpl}/out/{tpl}.phased.vcf.gz".format(tpl=self.prepare_tpl)

        if self.cfg.genes_of_interest_file:
            yield "genes", self.cfg.genes_of_interest_file
        if self.proteome_file:
            yield "peptides", self.proteome_file

    def _get_output_files_pileup(self):
        return {"vcf": "work/{tpl}/out/{tpl}.pileup.vcf.gz".format(tpl=self.prepare_tpl)}

    def _get_output_files_add_expression(self):
        return {"vcf": "work/{tpl}/out/{tpl}.with_expression.vcf.gz".format(tpl=self.prepare_tpl)}

    def _get_output_files_pvacseq(self):
        return self._get_output_files_run()

    def _get_args_pileup(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.pileup.model_dump(by_alias=True))

        return {
            "tumor_sample": self.get_sample_names(wildcards)["tumor_sample"],
            "extra_args": args["extra_args"],
        }

    def _get_args_add_expression(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.config.quantification.model_dump(by_alias=True))

        sample_names = self.get_sample_names(wildcards)
        return {
            "format": args["tool_gene_expression_quantification"],
            "extra_args": args["extra_args"],
            "tumor_sample": sample_names["tumor_sample"],
            "normal_sample": sample_names["normal_sample"],
        }

    def _get_args_pvacseq(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        n_threads = args["n_threads"]

        samples = self.get_sample_names(wildcards)

        extra_args = args["extra_args"]

        if self.cfg.ml_predictions:
            extra_args += [
                "--run-ml-predictions",
                f"--ml-threshold-accept {self.cfg.ml_predictions.accept}",
                f"--ml-threshold-reject {self.cfg.ml_predictions.reject}",
            ]

        excluded = ["container", "alleles"]
        for _, mhc_class_fn in self.SUBDIRECTORIES:
            for k, _ in self.FILE_EXTENSIONS:
                excluded.append(f"{k}.{mhc_class_fn}")

        return {
            "normal_sample": samples.get("normal_sample", None),
            "tumor_sample": samples["tumor_sample"],
            "algorithms": args["algorithms"],
            "lengths": {
                "class_i": args["class_i_epitope_length"],
                "class_ii": args["class_ii_epitope_length"],
            },
            "n_threads": n_threads,
            "exclude_bind": excluded,
            "extra_args": extra_args,
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

        input_files["alleles"] = f"work/{self.prepare_tpl}/out/{self.prepare_tpl}.hla_types.txt"

        if self.cfg.genes_of_interest_file:
            input_files["genes"] = self.cfg.genes_of_interest_file
        if self.proteome_file:
            input_files["peptides"] = self.proteome_file

        return input_files

    def _get_output_files_pvacfuse(self):
        return self._get_output_files_run()

    def _get_args_pvacfuse(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        n_threads = args["n_threads"]

        samples = self.get_sample_names(wildcards)

        excluded = ["container", "alleles"]
        for _, mhc_class_fn in self.SUBDIRECTORIES:
            for k, _ in self.FILE_EXTENSIONS:
                if k != "json":
                    excluded.append(f"{k}.{mhc_class_fn}")

        return {
            "tumor_sample": samples["tumor_sample"],
            "algorithms": args["algorithms"],
            "lengths": {
                "class_i": args["class_i_epitope_length"],
                "class_ii": args["class_ii_epitope_length"],
            },
            "n_threads": n_threads,
            "exclude_bind": excluded,
            "extra_args": args["extra_args"],
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
            time="03:59:59",
            memory="16G",
        ),
        "pvacsplice": ResourceUsage(
            threads=16,
            time="23:59:59",
            memory="64G",
        ),
        "workaround": ResourceUsage(
            threads=1,
            time="01:00:00",
            memory="4G",
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

        yield "alleles", f"work/{self.prepare_tpl}/out/{self.prepare_tpl}.hla_types.txt"

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

    def _get_args_pvacsplice(self, wildcards: Wildcards) -> dict[str, str]:
        args = dict(self.cfg.model_dump(by_alias=True))

        n_threads = args["n_threads"]

        samples = self.get_sample_names(wildcards)

        excluded = ["container", "alleles"]
        for _, mhc_class_fn in self.SUBDIRECTORIES:
            for k, _ in self.FILE_EXTENSIONS:
                if k != "json":
                    excluded.append(f"{k}.{mhc_class_fn}")

        return {
            "normal_sample": samples["normal_sample"],
            "tumor_sample": samples["tumor_sample"],
            "algorithms": args["algorithms"],
            "lengths": {
                "class_i": args["class_i_epitope_length"],
                "class_ii": args["class_ii_epitope_length"],
            },
            "n_threads": n_threads,
            "exclude_bind": excluded,
            "extra_args": args["extra_args"],
        }


class PhasingStepPart(BaseStepPart):
    """
    Phase somatic with germline variants using obsolete GATK, for pVACtools only.

    TODO: A better solution sould be developed using current GATK tools
    """

    name = "phasing"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, time="23:59:59", memory="32G")}

    def __init__(self, parent):
        super().__init__(parent)
        self.prepare_tpl = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            self.prepare_tpl += ".filtered"
        self.prepare_tpl += ".{tumor_dna}"

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
        tpl = "work/{tpl}/log/{action}.{{tumor_dna}}".format(tpl=self.prepare_tpl, action=self.name)
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


class PostProcessStepPart(GenericNeoepitopeStepPart):
    """base class to run netchop & netstab after pVACseq, pVACsplice or pVACfuse"""

    actions = ("run",)

    #: Resources
    default_resource_usage = ResourceUsage(threads=8, time="143:59:59", memory="32G")

    def __init__(self, parent):
        super().__init__(parent)
        prefix = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            postfix = "filtered.{tumor_dna}"
        else:
            postfix = "{tumor_dna}"
        self.prepare_tpl = f"{prefix}.{postfix}"
        self.input_tpl = f"{prefix}.{{tool}}.{postfix}"
        self.output_tpl = f"{prefix}.{{tool,pvacseq|pvacsplice|pvacfuse}}.{postfix}"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_log_files(self, tpl: str) -> dict[str, str]:
        """Return mapping of log files."""
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
        return self.default_resource_usage

    def _get_mhc_class(self, wildcards: Wildcards) -> MHC_CLASS:
        for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
            if wildcards.mhc_class_d == mhc_class.dirname:
                return mhc_class
        raise InvalidConfigurationException(f"Unknown MHC class {wildcards.mhc_class_d}")


class NetChopStepPart(PostProcessStepPart):
    """Run netchop"""

    name = "netchop"

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, str]:
        inputs = {"netchop": self.config.get(wildcards.tool).get("net_chop").get("path_netchop")}

        all_epitopes = (
            "all_epitopes"
            if self.config.get(wildcards.tool).get("net_chop").get("all_epitopes")
            else "filtered"
        )
        fn = os.path.join(
            "work",
            self.input_tpl,
            "out",
            "{mhc_class_d}",
            "{tumor_dna}.{mhc_class_fn}.{all_epitopes}.tsv",
        )
        inputs["epitopes"] = fn.format(all_epitopes=all_epitopes, **wildcards)
        # fn = os.path.join(
        #     "work",
        #     self.input_tpl,
        #     "out",
        #     "{mhc_class_d}",
        #     "{tumor_dna}.{mhc_class_fn}.fasta",
        # )
        # inputs["sequences"] = fn.format(**wildcards)

        return inputs

    def get_output_files(self, action: str) -> dict[str, Any]:
        self._validate_action(action)
        return {
            self.name: os.path.join(
                "work",
                self.output_tpl,
                "out",
                "{mhc_class_d}",
                f"{{tumor_dna}}.{{mhc_class_fn}}.{self.name}.tsv",
            )
        }

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        cfg: NetChopModel = self.config.get(wildcards.tool).get("net_chop")
        return {"tool": wildcards.tool, "method": cfg.method, "threshold": cfg.threshold}

    def get_log_file(self, action: str) -> dict[str, str]:
        self._validate_action(action)
        return self._get_log_files(
            f"work/{self.output_tpl}/log/{self.name}.{{mhc_class_d}}_{{mhc_class_fn}}.{{tumor_dna}}"
        )


class NetStabStepPart(PostProcessStepPart):
    """Run netMHCstabpan"""

    name = "netstab"
    actions = ("run",)

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, str]:
        if container := self.config.get(wildcards.tool).get("path_container"):
            pass
        else:
            container = "work/containers/out/pvactools.sif"

        if self.config.get(wildcards.tool).get("netmhc_stab").get("enabled"):
            ext = "netchop.tsv"
        else:
            ext = (
                "all_epitopes.tsv"
                if self.config.get(wildcards.tool).get("netmhc_stab").get("all_candidates")
                else "filtered.tsv"
            )
        fn = os.path.join(
            "work",
            self.input_tpl,
            "out",
            "MHC_Class_I",
            "{tumor_dna}.MHC_I." + ext,
        )

        return {
            "container": container,
            "epitopes": fn.format(**wildcards),
            "alleles": f"work/{self.prepare_tpl}/out/{self.prepare_tpl}.hla_types.txt",
        }

    def get_output_files(self, action: str) -> dict[str, Any]:
        self._validate_action(action)
        return {
            self.name: os.path.join(
                "work",
                self.output_tpl,
                "out",
                "MHC_Class_I",
                f"{{tumor_dna}}.MHC_I.{self.name}.tsv",
            )
        }

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "tool": wildcards.tool,
            "lengths": self.config.get(wildcards.tool).get("class_i_epitope_length"),
            "exclude_bind": ["container", "alleles"],
        }

    def get_log_file(self, action: str) -> dict[str, str]:
        self._validate_action(action)
        return self._get_log_files(f"work/{self.output_tpl}/log/{self.name}.{{tumor_dna}}")


class ProteomeStepPart(BaseStepPart):
    """Create personalised proteome"""

    name = "proteome"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, time="03:59:59", memory="24G")}

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        cfg: ProteomeModel = self.config.proteome
        if cfg.path_germline_variants:
            tpl = f"{cfg.tool_ngs_mapping}.{cfg.tool_germline_variant_calling}"
            if cfg.tool_germline_variant_annotation:
                tpl += f".{cfg.tool_germline_variant_annotation}"
            if cfg.is_filtered:
                tpl += ".filtered"
            tpl += ".{normal_dna}".format(normal_dna=self.parent.tumor_dna[wildcards.tumor_dna])
            germline_variant = self.parent.sub_workflows["germline_variant"]
            yield "vcf", germline_variant(f"output/{tpl}/out/{tpl}.vcf.gz")
        if cfg.external_proteome:
            yield "path_proteome", cfg.external_proteome
        yield "reference", self.w_config.static_data_config.reference.path
        yield "features", self.w_config.static_data_config.features.path

    def get_output_files(self, action: str) -> dict[str, str]:
        match action:
            case "run":
                if self.config.proteome.path_germline_variants:
                    tpl = "{mapper}.{caller}.{annotator}"
                    if self.config.is_filtered:
                        tpl += ".filtered"
                    tpl += ".{tumor_dna}"
                    return {"proteome": f"work/{tpl}/out/{tpl}.proteome.fa.gz"}
                else:
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
        if self.config.proteome.path_germline_variants:
            tpl = "{mapper}.{caller}.{annotator}"
            if self.config.is_filtered:
                tpl += ".filtered"
            tpl += ".{tumor_dna}"
            tpl = f"work/{tpl}/log/proteome.{{tumor_dna}}"
        else:
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


class HlaTypesStepPart(BaseStepPart):
    """Read & process HLA typing output"""

    name = "hla_types"
    actions = ("run",)

    #: Resources
    default_resource_usage = {"run": ResourceUsage(threads=1, time="00:59:59", memory="4G")}

    def __init__(self, parent):
        super().__init__(parent)
        self.path = "{mapper}.{caller}.{annotator}"
        if self.config.is_filtered:
            self.path += ".filtered"
        self.path += ".{tumor_dna}"

    def get_input_files(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def _get_input_files_run(self, wildcards: Wildcards) -> list[str]:
        normal_dna = self.parent.tumor_dna.get(wildcards.tumor_dna, None)
        tumor_rna = self.parent.tumor_rna.get(wildcards.tumor_dna, None)

        hla_typing = self.parent.sub_workflows.get("hla_typing")
        tpl = "output/{tool}.{library_name}/out/{tool}.{library_name}.json"
        hla_files = list()

        for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
            for tool in self.config.tools_hla_typing.get("dna", {}).get(mhc_class.name, []):
                prefix = self.w_config.step_config.get("hla_typing").get(tool).get("mapper", "")
                if prefix:
                    prefix += "."
                fn = tpl.format(tool=prefix + tool, library_name=wildcards.tumor_dna)
                hla_files.append(hla_typing(fn))
                if normal_dna:
                    fn = tpl.format(tool=prefix + tool, library_name=normal_dna)
                    hla_files.append(hla_typing(fn))
            if not tumor_rna:
                continue
            for tool in self.config.tools_hla_typing.get("rna", {}).get(mhc_class.name, []):
                prefix = self.w_config.step_config.get("hla_typing").get(tool).get("mapper", "")
                if prefix:
                    prefix += "."
                fn = tpl.format(tool=prefix + tool, library_name=tumor_rna)
                hla_files.append(hla_typing(fn))

        return sorted(hla_files)

    def get_output_files(self, action) -> dict[str, str]:
        self._validate_action(action)
        return {"hla_types": os.path.join("work", self.path, "out", self.path + ".hla_types.txt")}

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {}

    def get_log_file(self, action) -> dict[str, str]:
        self._validate_action(action)
        tpl = os.path.join("work", self.path, "log", "hla_types.{tumor_dna}")
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
        if cfg.get("proteome", {}).get("enabled", False) and cfg.get("proteome", {}).get(
            "path_germline_variants", None
        ):
            if cfg.get("proteome").get("tool_variant_annotation") == GermlineVariantStep.CALL:
                from snappy_pipeline.workflows.germline_variant_calling import (
                    GermlineVariantCallingWorkflow,
                )

                previous_steps.append(GermlineVariantCallingWorkflow)
            elif cfg.get("proteome").get("germline_variant_step") == GermlineVariantStep.FILTER:
                from snappy_pipeline.workflows.germline_variant_filtration import (
                    GermlineVariantFiltrationWorkflow,
                )

                previous_steps.append(GermlineVariantFiltrationWorkflow)
            else:
                from snappy_pipeline.workflows.germline_variant_annotation import (
                    GermlineVariantAnnotationWorkflow,
                )

                previous_steps.append(GermlineVariantAnnotationWorkflow)
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
                HlaTypesStepPart,
                PvacToolsStepPart,
                PvacSeqStepPart,
                PvacFuseStepPart,
                PvacSpliceStepPart,
                PhasingStepPart,
                NetChopStepPart,
                NetStabStepPart,
                ProteomeStepPart,
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
        if self.config.proteome.enabled and self.config.proteome.path_germline_variants:
            self.register_sub_workflow(
                self.config.proteome.germline_variant_step,
                self.config.proteome.path_germline_variants,
                "germline_variant",
            )
        if "pvacfuse" in self.config.tools:
            self.register_sub_workflow(
                "somatic_gene_fusion_calling",
                self.config.pvacfuse.path_somatic_gene_fusion_calling,
            )

    @listify
    def get_result_files(self):
        log_exts = ("log", "conda_list.txt", "conda_info.txt")
        hash_exts = ("", ".md5")
        exts = ["filtered.tsv", "all_epitopes.tsv", "all_epitopes.aggregated.tsv"]

        mappers = self.w_config.step_config["ngs_mapping"]["tools"]["dna"]
        callers = self.w_config.step_config["somatic_variant_calling"]["tools"]
        annotators = self.w_config.step_config["somatic_variant_annotation"]["tools"]

        library_prefix = "{mapper}.{caller}.{annotator}.{tool_name}"
        if self.config.is_filtered:
            library_prefix += ".filtered"

        tumor_samples = self.sample_table[
            (self.sample_table["extractionType"] == ExtractionType.DNA)
            & (self.sample_table["isTumor"])
        ]["ngs_library"]

        has_class_i = (
            len(
                self.config.get("tools_hla_typing").get("dna", {}).get("class_i", [])
                + self.config.get("tools_hla_typing").get("rna", {}).get("class_i", [])
            )
            > 0
        )
        has_class_ii = (
            len(
                self.config.get("tools_hla_typing").get("dna", {}).get("class_ii", [])
                + self.config.get("tools_hla_typing").get("rna", {}).get("class_ii", [])
            )
            > 0
        )

        for tool_name in self.config.tools:
            tool = self.sub_steps[tool_name]

            tool_has_class_i = has_class_i and any(
                map(
                    lambda a: a in ClassIAlgorithm,
                    self.config.get(tool_name).get("algorithms", []),
                )
            )
            tool_has_class_ii = has_class_ii and any(
                map(
                    lambda a: a in ClassIIAlgorithm,
                    self.config.get(tool_name).get("algorithms", []),
                )
            )
            tool_has_class_i = tool_has_class_i and self.config.get(tool_name).get(
                "class_i_epitope_length", []
            )
            tool_has_class_ii = tool_has_class_ii and self.config.get(tool_name).get(
                "class_ii_epitope_length", []
            )

            for tumor_dna in tumor_samples:
                if tool.require_rna and self.tumor_rna.get(tumor_dna, None) is None:
                    continue

                d = f"output/{library_prefix}.{tumor_dna}"
                fn = f"out/combined/{tumor_dna}.Combined.{{ext}}"
                yield from expand(
                    d + "/" + fn,
                    mapper=mappers,
                    caller=callers,
                    annotator=annotators,
                    tool_name=[tool_name],
                    ext=exts,
                )

                fn = f"log/{tool_name}.{tumor_dna}.{{log_ext}}{{hash_ext}}"
                yield from expand(
                    d + "/" + fn,
                    mapper=mappers,
                    caller=callers,
                    annotator=annotators,
                    tool_name=[tool_name],
                    log_ext=log_exts,
                    hash_ext=hash_exts,
                )

                if self.config.get(tool_name).get("net_chop").get("enabled"):
                    if tool_has_class_i:
                        fn = f"out/MHC_Class_I/{tumor_dna}.MHC_I.netchop.tsv"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                        )

                        fn = f"log/netchop.MHC_Class_I_MHC_I.{tumor_dna}.{{log_ext}}{{hash_ext}}"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                            log_ext=log_exts,
                            hash_ext=hash_exts,
                        )

                    if tool_has_class_ii:
                        fn = f"out/MHC_Class_II/{tumor_dna}.MHC_II.netchop.tsv"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                        )

                        fn = f"log/netchop.MHC_Class_II_MHC_II.{tumor_dna}.{{log_ext}}{{hash_ext}}"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                            log_ext=log_exts,
                            hash_ext=hash_exts,
                        )

                if self.config.get(tool_name).get("netmhc_stab").get("enabled"):
                    if tool_has_class_i:
                        fn = f"out/MHC_Class_I/{tumor_dna}.MHC_I.netstab.tsv"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                        )

                        fn = f"log/netstab.{tumor_dna}.{{log_ext}}{{hash_ext}}"
                        yield from expand(
                            d + "/" + fn,
                            mapper=mappers,
                            caller=callers,
                            annotator=annotators,
                            tool_name=[tool_name],
                            log_ext=log_exts,
                            hash_ext=hash_exts,
                        )

    def check_config(self):
        hla_typing_config = self.w_config.step_config.get("hla_typing", None)
        extraction_type = ExtractionType.DNA
        for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
            tool = self.config.tools_hla_typing.get(extraction_type, {}).get(mhc_class.name, None)
            if tool and tool not in hla_typing_config.tools.get(extraction_type, []):
                raise MissingConfiguration(f"hla_typing tool {tool} not configured")
        extraction_type = ExtractionType.RNA
        if (self.config.pileup.enabled or self.config.quantification.enabled) or (
            "pvacfuse" in self.config.tools or "pvacsplice" in self.config.tools
        ):
            for mhc_class in (MHC_CLASS_I, MHC_CLASS_II):
                tool = self.config.tools_hla_typing.get(extraction_type, {}).get(
                    mhc_class.name, None
                )
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

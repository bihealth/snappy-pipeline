# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_cnv_calling`` step

This step allows for the detection of CNV events for cancer samples from targeted sequenced (e.g.,
exomes or large panels) or whole genome sequencing. Panel sequencing is not implemented yet, it might be in a later release.
The wrapped tools start from the aligned reads (thus off ``ngs_mapping``) and generate CNV calls for somatic variants.

The wrapped tools implement different strategies.  Some work "reference free" and just use the
somatic BAM files for their input, some work in "matched cancer normal mode" and need the cancer
and normal BAM files, and finally others use a set of non-cancer BAM files for their background (the panel of normals).

Some tools may also use germline & somatic variants to estimate allele-specific copy number changes,
and resolve loss-of-heterozygocity. In this case, the small variants need to be computed separately from the ``somatic_variants_for_cnv`` step.

Finally, estimation of tumor purity and ploidy can sometimes improve considerably the reliability of CNV calls.
Some tools internally perform this estimation, possibly using germline & somatic small variants,
but others don't offer this possibility. In the latter case, these estimates might be provided by the
``somatic_purity_ploidy_estimate`` step, but in most cases, purity & ploidy values are so tightly connect to
the CNV calling procedure, that it is impossible to isolate this computation completely.
Therefore, purity tools have been included in the ``somatic_cnv_calling`` step, and their access restricted to
calling tools without internal purity estimation (such as ``cvnkit``).

==========
Step Input
==========

Gene somatic CNV calling starts off the aligned reads, i.e., ``ngs_mapping``.

Tools that use panel of normals can obtain their input in two different ways:

- A static file or files, from another cohort or from public datasets.
  In this case, the user is responsible to make sure that the data & methods used to create the panel are compatible to the cohort's.
- The ``panel_of_normals`` step.
  The panel will be created if necessary, using data from normal samples from the cohort.

When requested, the optional germline and somatic small variant calls are created in the ``somatic_variants_for_cnv`` step.
Once again, it is the responsability of the user to make sure that variants created in that way are suitable for CNV calling.

===========
Step Output
===========

TODO: The whole section of output needs revision. Main question is: what is the best format to encode CNAs?
``vcf`` is an possibility, the benefits are a (more or less) well-defined format, but the major drawback is
that (as far as I know), most CNV analysis tools (from ``R`` in particular) don't recognize this format (for CNAs).

Currently, there is no final formatting of the results into a format common to all tools.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_purity_ploidy_estimate.rst

============================
Available Somatic CNV Caller
============================

- Common to WES & WGS: ``ascat`` (later: ``cnvkit`` & ``sequenza``)
- WGS only: (later: ``Control-FREEC``)
- WES only: (later: ``PureCN``)

"""

from __future__ import annotations

import fnmatch

from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import InputFiles, Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.models.common import LibraryKitEntry, SexValue, SexOrigin
from snappy_pipeline.workflows.common.samplesheet import sample_sheets, tumor_to_normal_mapping
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.guess_sex import GuessSexWorkflow

from .model import SomaticCnvCalling as SomaticCnvCallingConfigModel
from .model import Ascat as AscatConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Tools for estimating purity and ploidy.
PURITY_PLOIDY_TOOLS = "ascat"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticCnvCallingConfigModel.default_config_yaml_string()


class SomaticCnvCallingStepPart(BaseStepPart):
    def _get_target_interval_mapping_path(self, library_name: str) -> str | None:
        if self.parent.libraryType == "WGS":
            return None
        libraryKit = self.parent.table.at[library_name, "libraryKit"]
        assert libraryKit in self.parent.kit_name_to_path, (
            f"No mapping for library kit '{libraryKit}'"
        )
        return self.parent.kit_name_to_path[libraryKit]

    def _get_guess_sex_file(self, wildcards: Wildcards) -> str:
        guess_sex = self.parent.sub_workflows["guess_sex"]
        base_path = (
            "output/{mapper}.{tool}.{library_name}/out/{mapper}.{tool}.{library_name}".format(
                mapper=wildcards.mapper,
                tool=self.config.sex.guess_sex_tool,
                library_name=wildcards.library_name,
            )
        )
        return guess_sex(base_path + ".txt")

    def _get_sex(self, wildcards: Wildcards | None = None) -> SexValue:
        match self.config.sex.source:
            case SexOrigin.AUTOMATIC:
                sex = self._read_sex(self._get_guess_sex_file(wildcards))
            case SexOrigin.SAMPLESHEET:
                sex = SexValue(
                    self.parent.table.at[wildcards.library_name, self.config.sex.column_name]
                )
            case SexOrigin.CONFIG:
                sex = self.config.sex.cohort
        return sex

    def _read_sex(self, filename: str) -> SexValue:
        sex: SexValue | None = None
        with open(filename, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if line in SexValue._value2member_map_:
                    if sex is not None and sex != SexValue(line):
                        raise ValueError(f"Multiple incompatible sex values in file '{filename}'")
                    sex = SexValue(line)
        assert sex is not None, f"No valid sex value in {filename}"
        return sex


class AscatStepPart(SomaticCnvCallingStepPart):
    """Somatic cnv calling using ASCAT."""

    # Name of the step.
    name = "ascat"

    # The actions for generating BAF/CNV files for the tumor and/or normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "alleleCounter",
        "prepareHTS",
        "run",
    )

    resource_usage = {
        "alleleCounter": ResourceUsage(threads=1, time="4:00:00", memory=f"{4 * 1024}M"),
        "prepareHTS": ResourceUsage(threads=1, time="4:00:00", memory=f"{16 * 1024}M"),
        "run": ResourceUsage(threads=8, time="2-00:00:00", memory=f"{10 * 1024 * 8}M"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: AscatConfigModel = self.config.get(self.name)
        self.chromosome_names = AscatStepPart._get_chromosome_names(
            ignored=set(getattr(self.config, "ignore_chroms", []))
            | set(getattr(self.cfg, "ignore_chroms", []))
        )

    @listify
    @staticmethod
    def _get_chromosome_names(
        all_chromosomes: list[str] = list(map(str, range(1, 23))) + ["X"], ignored: set[str] = {}
    ):
        """
        Internally, ASCAT requires chromosome names as 1, ..., 22 & X (CNV on Y are not called).
        These names (without prefix) are used when calling ASCAT R functions, and in
        filenames read and/or produced by ASCAT.
        Note that the ASCAT files are ..._chr<1,2,...,X>.txt, regardless whether the actual
        chromosome name has a chr prefix or not.

        This imposes some restrictions to the set of ignored chromosomes:
        if the user wants  to exclude some chromosome, she must use the un-prefixed name in
        the ``ignore_chroms`` configuration, even if the data is mapped against a prefixed
        genome. We suggest that in this case, the user prefers the tool's ``ignore_chroms`` option,
        rather than the step's.
        """
        for contig_name in all_chromosomes:
            keep = True
            for pattern in ignored:
                if fnmatch.fnmatch(contig_name, pattern):
                    keep = False
                    break
            if keep:
                yield contig_name

    def _get_gender(self, wildcards: Wildcards) -> str:
        return "XX" if self._get_sex(wildcards) == SexValue.FEMALE else "XY"

    def get_input_files(self, action: str):
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def get_output_files(self, action: str):
        self._validate_action(action)
        if action == "alleleCounter":
            tpl = "work/{mapper}.ascat.{tumor_name}/tmp/{library_name}_alleleFrequencies_chr{chrom_name}."
            yield "freq", tpl + "txt"
            yield "vcf", tpl + "vcf.gz"
        else:
            tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
            if action == "prepareHTS":
                output_files = {
                    "tumor_logr": tpl + "Tumor_LogR.txt",
                    "tumor_baf": tpl + "Tumor_BAF.txt",
                    "normal_logr": tpl + "Germline_LogR.txt",
                    "normal_baf": tpl + "Germline_BAF.txt",
                }
            elif action == "run":
                output_files = {
                    "RData": tpl + "RData",
                    "goodness_of_fit": tpl + "goodness_of_fit.txt",
                    "na": tpl + "na.txt",
                    "nb": tpl + "nb.txt",
                    "ploidy": tpl + "ploidy.txt",
                    "purity": tpl + "purity.txt",
                    "segments": tpl + "segments.txt",
                    "segments_raw": tpl + "segments_raw.txt",
                }
            for k, v in output_files.items():
                yield k, v
                yield k + "_md5", v + ".md5"

    def get_result_files(self) -> list[str]:
        tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
        results = [
            tpl + "goodness_of_fit.txt",
            tpl + "na.txt",
            tpl + "nb.txt",
            tpl + "ploidy.txt",
            tpl + "purity.txt",
            tpl + "segments.txt",
        ]
        tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}."
        for action in ("prepareHTS", "run"):
            for ext in ("conda_info.txt", "conda_list.txt", "log", "R"):
                results.append(tpl + f"{action}.{ext}")
        return results

    @dictify
    def get_log_file(self, action: str):
        self._validate_action(action)
        if action == "alleleCounter":
            tpl = (
                "work/{mapper}.ascat.{tumor_name}/log/alleleCounter.{library_name}.chr_{chrom_name}"
            )
            yield "script", tpl + ".alleleCounter.sh"
            yield "script_md5", tpl + ".alleleCounter.sh.md5"
        else:
            tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}"
            k = "script"
            v = tpl + "." + action + ".R"
            yield k, v
            yield k + "_md5", v + ".md5"
        for ext in ("conda_info.txt", "conda_list.txt", "log"):
            k = ext.replace(".txt", "")
            v = tpl + "." + action + "." + ext
            yield k, v
            yield k + "_md5", v + ".md5"

    def get_args(self, action: str):
        self._validate_action(action)
        return getattr(self, "_get_args_{}".format(action))

    def _get_input_files_alleleCounter(self, wildcards: Wildcards) -> dict[str, str]:
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "reference": self.w_config.static_data_config.reference.path,
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
            "loci": self.cfg.allele_counter.loci_prefix + f"{wildcards.chrom_name}.txt",
        }
        return input_files

    def _get_args_alleleCounter(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        args = {
            "seed": self.cfg.seed,
            "max-depth": self.cfg.allele_counter.max_depth,
            "min-MQ": self.cfg.allele_counter.min_MQ,
            "min-BQ": self.cfg.allele_counter.min_BQ,
            "skip-any-set": self.cfg.allele_counter.skip_any_set.value,
            "skip-any-unset": self.cfg.allele_counter.skip_any_unset.value,
        }

        return args

    def _get_input_files_prepareHTS(self, wildcards: Wildcards) -> dict[str, str]:
        input_files = {
            "reference": self.w_config.static_data_config.reference.path,
            "tumor_alleleFrequencies": [
                "work/{mapper}.ascat.{tumor_name}/tmp/{library_name}_alleleFrequencies_chr{chrom_name}.txt".format(
                    mapper=wildcards.mapper,
                    tumor_name=wildcards.library_name,
                    library_name=wildcards.library_name,
                    chrom_name=chrom_name,
                )
                for chrom_name in self.chromosome_names
            ],
            "normal_alleleFrequencies": [
                "work/{mapper}.ascat.{tumor_name}/tmp/{library_name}_alleleFrequencies_chr{chrom_name}.txt".format(
                    mapper=wildcards.mapper,
                    tumor_name=wildcards.library_name,
                    library_name=self.parent.normals[wildcards.library_name],
                    chrom_name=chrom_name,
                )
                for chrom_name in self.chromosome_names
            ],
        }

        tpl = self.cfg.allele_counter.allele_prefix + "{chrom_name}.txt"
        input_files["alleles"] = [
            tpl.format(chrom_name=chrom_name, **wildcards) for chrom_name in self.chromosome_names
        ]

        if self.parent.libraryType != "WGS":
            input_files["baits"] = self._get_target_interval_mapping_path(wildcards.library_name)

        if self.config.sex.source == SexOrigin.AUTOMATIC:
            input_files["sex"] = self._get_guess_sex_file(wildcards)

        if f := self.cfg.allele_counter.path_probloci_file:
            input_files["probloci_file"] = f

        return input_files

    def _get_args_prepareHTS(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {
            "genomeVersion": self.cfg.genomeVersion,
            "minCounts": self.cfg.allele_counter.minCounts,
            "seed": self.cfg.seed,
            "chrom_names": self.chromosome_names,
            "gender": self._get_gender(wildcards),
            "tumorname": wildcards.library_name,
            "normalname": self.parent.normals[wildcards.library_name],
        }

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, str]:
        """Return input files for actually running ASCAT."""
        input_files = {"GCcontent": self.cfg.path_gc_content, "reptiming": self.cfg.path_reptiming}

        tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.".format(
            **wildcards
        )

        input_files["tumor_logr"] = tpl + "Tumor_LogR.txt"
        input_files["tumor_baf"] = tpl + "Tumor_BAF.txt"
        input_files["normal_logr"] = tpl + "Germline_LogR.txt"
        input_files["normal_baf"] = tpl + "Germline_BAF.txt"

        if self.config.sex.source == SexOrigin.AUTOMATIC:
            input_files["sex"] = self._get_guess_sex_file(wildcards)

        return input_files

    def _get_args_run(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {
            "genomeVersion": self.cfg.genomeVersion,
            "seed": self.cfg.seed,
            "gender": self._get_gender(wildcards),
            "y_limit": self.cfg.y_limit,
            "penalty": self.cfg.advanced.penalty,
            "gamma": self.cfg.advanced.gamma,
            "min_purity": self.cfg.advanced.min_purity,
            "max_purity": self.cfg.advanced.max_purity,
            "min_ploidy": self.cfg.advanced.min_ploidy,
            "max_ploidy": self.cfg.advanced.max_ploidy,
            "rho_manual": self.cfg.advanced.rho_manual,
            "psi_manual": self.cfg.advanced.psi_manual,
            "chrom_names": self.chromosome_names,
        }


class SomaticCnvCallingWorkflow(BaseStep):
    """Perform somatic CNV calling"""

    #: Workflow name
    name = "somatic_cnv_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticCnvCallingConfigModel,
            previous_steps=(NgsMappingWorkflow, GuessSexWorkflow),
        )

        self.table = sample_sheets(self.sheets)

        for column_name in ("isTumor", "extractionType", "libraryType"):
            assert column_name in self.table.columns, (
                f"Mandatory column '{column_name}' is missing from samplesheet"
            )

        self.table = self.table[self.table["extractionType"] == "DNA"]

        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        if self.config.sex.source == SexOrigin.AUTOMATIC:
            self.register_sub_workflow("guess_sex", self.config.sex.path_guess_sex)

        self.register_sub_step_classes((AscatStepPart, LinkOutStepPart))

        self.normals = tumor_to_normal_mapping(self.table)
        self.libraryType = self._get_library_type()
        if self.libraryType != "WGS":
            self.kit_name_to_path = self._get_defined_library_kits(
                self.config.path_target_interval_list_mapping
            )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        if self.table.shape[0] == 0:
            self.logger.warning("No sample to process")
        for sub_step in self.sub_steps.values():
            tool = sub_step.name
            tool_config = getattr(self.config, tool, None)
            if tool_config is None or self.table.shape[0] == 0:
                continue
            for ngs_library in self.table[self.table["isTumor"]].index.to_list():
                files_in_work = sub_step.get_result_files()
                for file_in_work in files_in_work:
                    file_in_output = file_in_work.format(
                        mapper=self.config.tool_ngs_mapping, library_name=ngs_library
                    ).replace("work/", "output/", 1)
                    yield file_in_output
                    yield file_in_output + ".md5"

    # TODO: generalize & put in common.samplesheet
    def _get_library_type(self) -> str:
        is_wes = 0
        is_wgs = 0
        is_other = 0
        for ngs_library in self.table.index.to_list():
            libraryType = self.table.at[ngs_library, "libraryType"]
            if libraryType == "WGS":
                is_wgs = 1
            elif libraryType == "WES":
                is_wes = 1
            else:
                is_other = 1
        if is_wes + is_wgs + is_other > 1:
            raise ValueError(
                "The current implementation doesn't support multiple library types in the same cohort"
            )
        if is_other == 1:
            raise ValueError(f"Unimplemented library type '{libraryType}'")
        return "WES" if is_wes else "WGS"

    def _get_defined_library_kits(
        self, baits_definition: list[LibraryKitEntry] = []
    ) -> dict[str, str]:
        assert "libraryKit" in self.table.columns, (
            "Mandatory column 'libraryKit' is missing from samplesheet"
        )
        assert len(baits_definition) > 0, "Missing baits definition"
        return {bait.name: bait.path for bait in baits_definition}

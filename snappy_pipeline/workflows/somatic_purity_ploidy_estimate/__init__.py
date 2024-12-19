# -*- coding: utf-8 -*-
"""Implementation of purity and ploidy checking for somatic NGS samples

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_purity_ploidy_estimate.rst

"""

from collections import OrderedDict
from typing import Any

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, KEY_IS_TUMOR
from snakemake.io import InputFiles, Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.models.common import SexValue, SexOrigin
from snappy_pipeline.workflows.common.samplesheet import sample_sheets
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.guess_sex import GuessSexWorkflow

from .model import SomaticPurityPloidyEstimate as SomaticPurityPloidyEstimateConfigModel
from .model import Ascat as AscatConfigModel

__author__ = "Eric Blanc <eric.blanc@bih-charite.de>"

#: Tools for estimating purity and ploidy.
PURITY_PLOIDY_TOOLS = "ascat"

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticPurityPloidyEstimateConfigModel.default_config_yaml_string()


class AscatStepPart(BaseStepPart):
    """Estimation of purity and ploidy using ASCAT.

    Notes
    -----

    - Although we name the virtual probes "SNP${num}", they are not guaranteed
      to
    """

    # Name of the step.
    name = "ascat"

    # The actions for generating BAF/CNV files for the tumor and/nor normal
    # sample, and finally to run the ASCAT pipeline.
    actions = (
        "build_baf",
        "prepare_hts",
        "run",
    )

    resource_usage = {
        "build_baf": ResourceUsage(threads=2, time="4:00:00", memory=f"{4 * 1024}M"),
        "prepare_hts": ResourceUsage(threads=2, time="4:00:00", memory=f"{4 * 1024}M"),
        "run": ResourceUsage(threads=8, time="2-00:00:00", memory=f"{10 * 1024 * 8}M"),
    }

    # X should be present even for male cohorts (https://github.com/VanLoo-lab/ascat/issues/125)
    chromosome_names = list(map(str, range(1, 23))) + ["X"]

    def __init__(self, parent):
        super().__init__(parent)
        self.cfg: AscatConfigModel = self.config.get(self.name)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

        self.baits = {}
        if "libraryType" in self.parent.table.columns and "libraryKit" in self.parent.table.columns:
            for row in self.parent.table.itertuples():
                ngs_library = row.Index
                if row.libraryType == "WES":
                    libraryKit = row.libraryKit
                    found = False
                    for entry in self.cfg.path_target_interval_list_mapping:
                        if entry.name == libraryKit:
                            self.baits[ngs_library] = entry.path
                            found = True
                            break
                    assert (
                        found
                    ), f"No exome kit description matching {libraryKit} for library {ngs_library}"

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.library_name]
        return pair.normal_sample.dna_ngs_library.name

    def _get_guess_sex_file(self, wildcards: Wildcards) -> str:
        guess_sex = self.parent.sub_workflows["guess_sex"]
        base_path = (
            "output/{mapper}.{tool}.{library_name}/out/{mapper}.{tool}.{library_name}".format(
                mapper=wildcards["mapper"],
                tool=self.cfg.sex.guess_sex_tool,
                library_name=wildcards["library_name"],
            )
        )
        return guess_sex(base_path + ".txt")

    def _get_gender(self, library_name: str, wildcards: Wildcards | None = None) -> str:
        match self.cfg.sex.source:
            case SexOrigin.AUTOMATIC:
                sex = self._read_sex(self._get_guess_sex_file(wildcards))
            case SexOrigin.SAMPLESHEET:
                sex = self.parent.table.loc[library_name][self.cfg.sex.column_name]
            case SexOrigin.CONFIG:
                sex = self.cfg.sex.cohort
        if sex == SexValue.MALE:
            gender = "XY"
        elif sex == SexValue.FEMALE:
            gender = "XX"
        else:
            raise ValueError(f"Can't obtain sex for library {library_name}")
        return gender

    def get_input_files(self, action: str):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @dictify
    def get_output_files(self, action: str):
        # Validate action
        self._validate_action(action)
        tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}."
        match action:
            case "build_baf":
                tpl = "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}."
                output_files = {"vcf": tpl + "vcf.gz", "txt": tpl + "txt"}
            case "prepare_hts":
                output_files = {
                    "tumor_logr": tpl + "Tumor_LogR.txt",
                    "tumor_baf": tpl + "Tumor_BAF.txt",
                    "normal_logr": tpl + "Germline_LogR.txt",
                    "normal_baf": tpl + "Germline_BAF.txt",
                }
            case "run":
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
        for action in ("prepare_hts", "run"):
            tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}."
            for ext in ("conda_info.txt", "conda_list.txt", "log", "R"):
                results.append(tpl + f"{action}.{ext}")
        return results

    @dictify
    def get_log_file(self, action: str):
        # Validate action
        self._validate_action(action)
        if action == "build_baf":
            tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}.chr{chrom_name}."
        else:
            tpl = "work/{mapper}.ascat.{library_name}/log/{mapper}.ascat.{library_name}."
        for ext in ("conda_info.txt", "conda_list.txt", "log"):
            k = ext.replace(".txt", "")
            v = tpl + action + "." + ext
            yield k, v
            yield k + "_md5", v + ".md5"
        if action in ("build_baf",):
            ext = "sh"
        else:
            ext = "R"
        k = "script"
        v = tpl + action + "." + ext
        yield k, v
        yield k + "_md5", v + ".md5"

    def get_args(self, action: str):
        # Validate action
        self._validate_action(action)
        return getattr(self, "_get_args_{}".format(action))

    def _get_input_files_build_baf(self, wildcards: Wildcards) -> dict[str, str]:
        """Return input files for generating BAF file for either tumor or normal."""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        base_path = "output/{mapper}.{library_name}/out/{mapper}.{library_name}".format(**wildcards)
        input_files = {
            "reference": self.w_config.static_data_config.reference.path,
            "locii": self.cfg.allele_counter.loci_prefix + "{chrom_name}.txt".format(**wildcards),
            "bam": ngs_mapping(base_path + ".bam"),
            "bai": ngs_mapping(base_path + ".bam.bai"),
        }
        if wildcards["library_name"] in self.baits:
            input_files["baits"] = self.baits[wildcards["library_name"]]
        return input_files

    def _get_args_build_baf(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        return {
            "minCounts": self.cfg.allele_counter.minCounts,
            "min_map_qual": self.cfg.allele_counter.min_map_qual,
            "min_base_qual": self.cfg.allele_counter.min_base_qual,
            "max_coverage": self.cfg.allele_counter.max_coverage,
            "exclude_flags": self.cfg.allele_counter.exclude_flags,
            "include_flags": self.cfg.allele_counter.include_flags,
        }

    def _get_input_files_prepare_hts(self, wildcards: Wildcards) -> dict[str, str]:
        tpl = "work/{mapper}.ascat.{library_name}/out/AlleleCounts/{library_name}_alleleFrequencies_chr{chrom_name}.txt"
        input_files = {
            "tumorAlleleCounts": [
                tpl.format(
                    mapper=wildcards["mapper"],
                    library_name=wildcards["library_name"],
                    chrom_name=chrom_name,
                )
                for chrom_name in AscatStepPart.chromosome_names
            ],
        }

        if self.tumor_ngs_library_to_sample_pair.get(wildcards["library_name"], None) is not None:
            normal_library = self.get_normal_lib_name(wildcards)
            input_files["normalAlleleCounts"] = [
                tpl.format(
                    mapper=wildcards["mapper"],
                    library_name=normal_library,
                    chrom_name=chrom_name,
                )
                for chrom_name in AscatStepPart.chromosome_names
            ]

        tpl = self.cfg.allele_counter.allele_prefix + "{chrom_name}.txt"
        input_files["alleles"] = [
            tpl.format(chrom_name=chrom_name) for chrom_name in AscatStepPart.chromosome_names
        ]

        if self.cfg.sex.source == SexOrigin.AUTOMATIC:
            input_files["sex"] = self._get_guess_sex_file(wildcards)

        if self.cfg.allele_counter.path_probloci_file:
            input_files["probloci_file"] = self.cfg.allele_counter.path_probloci_file

        return input_files

    def _get_args_prepare_hts(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        args = {
            "genomeVersion": self.cfg.genomeVersion,
            "min_base_qual": self.cfg.allele_counter.min_base_qual,
            "min_map_qual": self.cfg.allele_counter.min_map_qual,
            "minCounts": self.cfg.allele_counter.minCounts,
            "seed": self.cfg.seed,
            "chrom_names": AscatStepPart.chromosome_names,
            "gender": self._get_gender(wildcards["library_name"], wildcards),
            "tumorname": wildcards["library_name"],
        }

        if self.tumor_ngs_library_to_sample_pair.get(wildcards["library_name"], None) is not None:
            args["normalname"] = self.get_normal_lib_name(wildcards)

        return args

    def _get_input_files_run(self, wildcards: Wildcards) -> dict[str, str]:
        """Return input files for actually running ASCAT."""
        input_files = {"GCcontent": self.cfg.path_gc_content, "reptiming": self.cfg.path_reptiming}

        tpl = "work/{mapper}.ascat.{library_name}/out/{mapper}.ascat.{library_name}.".format(
            **wildcards
        )

        input_files["tumor_logr"] = tpl + "Tumor_LogR.txt"
        input_files["tumor_baf"] = tpl + "Tumor_BAF.txt"
        if self.tumor_ngs_library_to_sample_pair.get(wildcards["library_name"], None) is not None:
            input_files["normal_logr"] = tpl + "Germline_LogR.txt"
            input_files["normal_baf"] = tpl + "Germline_BAF.txt"

        if self.cfg.sex.source == SexOrigin.AUTOMATIC:
            input_files["sex"] = self._get_guess_sex_file(wildcards)

        return input_files

    def _get_args_run(self, wildcards: Wildcards, input: InputFiles) -> dict[str, Any]:
        args = {
            "genomeVersion": self.cfg.genomeVersion,
            "seed": self.cfg.seed,
            "gender": self._get_gender(wildcards["library_name"], wildcards),
            "y_limit": self.cfg.y_limit,
            "penalty": self.cfg.advanced.penalty,
            "gamma": self.cfg.advanced.gamma,
            "min_purity": self.cfg.advanced.min_purity,
            "max_purity": self.cfg.advanced.max_purity,
            "min_ploidy": self.cfg.advanced.min_ploidy,
            "max_ploidy": self.cfg.advanced.max_ploidy,
            "rho_manual": self.cfg.advanced.rho_manual,
            "psi_manual": self.cfg.advanced.psi_manual,
            "chrom_names": AscatStepPart.chromosome_names,
        }
        if self.tumor_ngs_library_to_sample_pair.get(wildcards["library_name"], None) is None:
            args["platform"] = self.cfg.platform
        return args

    def _read_sex(self, filename: str) -> SexValue:
        sex: SexValue | None = None
        with open(filename, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#"):
                    continue
                if line in SexValue._value2member_map_:
                    sex = line
                    break
        assert sex is not None, f"No valid sex value in {filename}"
        return sex


class SomaticPurityPloidyEstimateWorkflow(BaseStep):
    """Perform purity and ploidy estimation"""

    #: Workflow name
    name = "somatic_purity_ploidy_estimate"

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
            config_model_class=SomaticPurityPloidyEstimateConfigModel,
            previous_steps=(NgsMappingWorkflow, GuessSexWorkflow),
        )

        self.table = sample_sheets(self.sheets)
        assert (
            KEY_IS_TUMOR in self.table.columns
        ), f"Mandatory column '{KEY_IS_TUMOR}' is missing from samplesheet"
        if "extractionType" in self.table.columns:
            self.table = self.table[self.table["extractionType"] == "DNA"]

        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        path_guess_sex = None
        for tool in self.config.tools:
            tool_config = getattr(self.config, tool)  # config should be present, already validated
            if getattr(tool_config, "sex", None) and tool_config.sex.source == SexOrigin.AUTOMATIC:
                if path_guess_sex is None:
                    path_guess_sex = tool_config.sex.path_guess_sex
                assert (
                    path_guess_sex == tool_config.sex.path_guess_sex
                ), "Multiple paths to the output of 'guess_sex' step"
        if path_guess_sex is not None:
            self.register_sub_workflow("guess_sex", path_guess_sex)

        self.register_sub_step_classes((AscatStepPart, LinkOutStepPart))

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
                if not getattr(row, KEY_IS_TUMOR) or row.libraryType not in ("WES", "WGS"):
                    continue
                files_in_work = sub_step.get_result_files()
                for file_in_work in files_in_work:
                    file_in_output = file_in_work.format(
                        mapper=self.config.tool_ngs_mapping, library_name=row.Index
                    ).replace("work/", "output/", 1)
                    yield file_in_output
                    yield file_in_output + ".md5"

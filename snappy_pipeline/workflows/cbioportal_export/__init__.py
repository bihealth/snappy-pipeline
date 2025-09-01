# -*- coding: utf-8 -*-
"""Implementation of the ``cbioportal_export`` step

This step takes as input a wide range of files previously produced by the
snappy_pipeline. It does the necessary data transformations to comply with the
formats expected by cBioPortal and writes out the necessary metadata and data
files.
"""

import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage

from .model import CbioportalExport as CbioportalExportConfigModel

# cbioportal meta data files
META_FILES = {
    "always_present": ["meta_study.txt", "meta_clinical_patient.txt", "meta_clinical_sample.txt"],
    "sequenced": ["meta_mutation_extended.txt"],
    "cna": ["meta_cna_gistic.txt", "meta_cna_log2.txt"],
    "rna_seq_mrna": ["meta_expression.txt"],
    "segment": ["meta_segment.txt"],
}

# files holding clinical data
CLINICAL_DATA_FILES = ["data_clinical_patients.txt", "data_clinical_samples.txt"]

# files with CNA results
CNA_DATA_FILES = ["data_cna_gistic.txt", "data_cna_log2.txt", "data_segment.txt"]

# case lists files
CASE_LIST_FILES = {
    "sequenced": {
        "filename": "all_cases_with_mutation_data.txt",
        "name": "Sequenced tumors",
        "description": "Tumors with somatic variant calls",
        "stable_id": "sequenced",
        "category": "all_cases_with_mutation_data",
    },
    "cna": {
        "filename": "all_cases_with_cna_data.txt",
        "name": "Tumors with CNA data",
        "description": "Tumors with somatic Copy Number Alteration calls",
        "stable_id": "cna",
        "category": "all_cases_with_cna_data",
    },
    "rna_seq_mrna": {
        "filename": "all_cases_with_mrna_rnaseq_data.txt",
        "name": "Tumors with expression data",
        "description": "Tumors with mRNA seq expression data",
        "stable_id": "rna_seq_mrna",
        "category": "all_cases_with_mrna_rnaseq_data",
    },
    "cnaseq": {
        "filename": "all_cases_with_mutation_and_cna_data.txt",
        "name": "Sequenced tumors with CNA",
        "description": "Tumors with somatic variant & CNA calls",
        "stable_id": "cnaseq",
        "category": "all_cases_with_mutation_and_cna_data",
    },
    "3way_complete": {
        "filename": "all_cases_with_mutation_and_cna_and_mrna_data.txt",
        "name": "Sequenced tumors with CNA & expression",
        "description": "Tumors with somatic variant calls, CNA calls & expression data",
        "stable_id": "3way_complete",
        "category": "all_cases_with_mutation_and_cna_and_mrna_data",
    },
}


DEFAULT_CONFIG = CbioportalExportConfigModel.default_config_yaml_string()


# ================================================================================================
#
# Abstract classes: one generic, and the other for substeps dealing with a single sample at a time
#
# ================================================================================================


class cbioportalExportStepPart(BaseStepPart):
    """Base class for all cBioPortal export step parts

    - Builds a mapping between the tumor library id and the normal/tumor sample pair
    - Provides a list of input files based on a template, and on the iteration over the sheets.
      The template can be changed by derived classes.
    """

    #: Class available actions
    actions = None

    #: Extraction type, must be instantiated in sub-class
    extraction_type = None

    #: input files template
    input_tpl = None

    #: Output_file name
    output_file = None

    def _yield_libraries(self):
        test_msg = 'WARNING: multiple test samples for sample "{}", test sample "{}" is ignored'
        lib_msg = 'WARNING: multiple libraries for sample "{}", library "{}" is ignored'
        samples = []
        for sheet in filter(is_not_background, self.parent.sheets):
            for donor in sheet.bio_entities.values():
                for biosample in donor.bio_samples.values():
                    if biosample.extra_infos["isTumor"]:
                        sample_name = biosample.name
                        for test_sample in biosample.test_samples.values():
                            if (
                                self.extraction_type
                                and test_sample.extra_infos["extractionType"]
                                != self.extraction_type
                            ):
                                continue
                            if sample_name in samples:
                                print(
                                    test_msg.format(sample_name, test_sample.name), file=sys.stderr
                                )
                                continue
                            samples.append(sample_name)
                            first = True
                            for lib in test_sample.ngs_libraries.values():
                                # pick only one library
                                if not first:
                                    print(lib_msg.format(sample_name, lib.name), file=sys.stderr)
                                    break
                                yield lib
                                first = False

    @dictify
    def get_input_files(self, action):
        """Return path of input files for merging"""
        # Validate action
        self._validate_action(action)
        assert self.input_tpl is not None
        for lib in self._yield_libraries():
            yield lib.test_sample.bio_sample.name, self.input_tpl.format(library_name=lib.name)

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        assert self.output_file is not None
        return self.output_file

    def get_log_file(self, action):
        """Return path to log files for all data files"""
        # Validate action
        self._validate_action(action)

        tpl = "work/log/{name}.{action}{ext}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        log_files = {}
        for key, ext in key_ext:
            log_files[key] = tpl.format(name=self.name, action=action, ext=ext)
            log_files[key + "_md5"] = log_files[key] + ".md5"
        return log_files

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        return dict(self.config)


class cbioportalVcf2MafStepPart(BaseStepPart):
    """Helper class for VCF2MAF step"""

    #: Step name
    name = "cbioportal_vcf2maf"

    #: Actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.name_pattern = None
        if self.config.is_filtered:
            self.name_pattern = "{mapper}.{caller}.{annotator}.filtered.{tumor_library}"
        else:
            self.name_pattern = "{mapper}.{caller}.{annotator}.{tumor_library}"
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    @dictify
    def get_output_files(self, action):
        """Return maf output file"""
        # Validate action
        self._validate_action(action)
        yield "maf", os.path.join("work/maf", self.name_pattern, "out", self.name_pattern + ".maf")

    @dictify
    def get_input_files(self, action):
        """Return input vcf for each output maf"""
        # Validate action
        self._validate_action(action)
        somatic_variant = self.parent.sub_workflows["somatic_variant"]
        tpl = somatic_variant(
            os.path.join("output", self.name_pattern, "out", self.name_pattern + ".vcf.gz")
        )
        yield "vcf", tpl

    @dictify
    def get_log_file(self, action):
        """Return path to log files for all data files"""
        # Validate action
        self._validate_action(action)

        tpl = os.path.join("work/maf/", self.name_pattern, "log", self.name_pattern)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        def args_function(wildcards):
            result = {
                "tumor_sample": wildcards.tumor_library,
                "normal_sample": self._get_normal_lib_name(wildcards),
                "tumor_id": self._get_tumor_bio_sample(wildcards),
                "normal_id": self._get_normal_bio_sample(wildcards),
                "somatic_variant_annotation_tool": self.config.somatic_variant_annotation_tool,
                "ncbi_build": self.config.vcf2maf.ncbi_build,
                "Center": self.config.vcf2maf.Center,
            }
            return result

        return args_function

    def _get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def _get_tumor_bio_sample(self, wildcards):
        """Return bio sample to tumor dna ngs library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.tumor_sample.dna_ngs_library.test_sample.bio_sample.name

    def _get_normal_bio_sample(self, wildcards):
        """Return normal bio sample to tumor dna ngs library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.test_sample.bio_sample.name

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="5120M",
        )


class cbioportalMutationsStepPart(cbioportalExportStepPart):
    """Merges sample MAF files into data_mutation_extended.txt"""

    #: Step name
    name = "cbioportal_mutations"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_mutation_extended.txt"

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = "{mapper}.{caller}.{annotator}."
        if self.config.is_filtered:
            name_pattern += "filtered.{{library_name}}"
        else:
            name_pattern += "{{library_name}}"
        tpl = os.path.join("work/maf", name_pattern, "out", name_pattern + "{ext}")
        self.input_tpl = tpl.format(
            mapper=self.config.mapping_tool,
            caller=self.config.somatic_variant_calling_tool,
            annotator=self.config.somatic_variant_annotation_tool,
            ext=".maf",
        )


class cbioportalCns2CnaStepPart(BaseStepPart):
    """Helper class to extract gene-based log2 & copy number call from a cns file"""

    #: Step name
    name = "cbioportal_cns2cna"

    #: Actions
    actions = ("run",)

    @dictify
    def get_input_files(self, action):
        """Return the library"""
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{tumor_library}"
        yield "features", self.parent.w_config.static_data_config.features.path
        yield (
            "DNAcopy",
            os.path.join(
                self.config.copy_number_alteration.path_copy_number,
                "output",
                name_pattern,
                "out",
                name_pattern + "_dnacopy.seg",
            ),
        )

    @dictify
    def get_output_files(self, action):
        """Return maf output file"""
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{tumor_library}"
        yield "cna", os.path.join("work/cna", name_pattern, "out", name_pattern + ".cna")

    @dictify
    def get_log_file(self, action):
        """Return path to log files for all data files"""
        # Validate action
        self._validate_action(action)
        name_pattern = "{mapper}.{caller}.{tumor_library}"
        tpl = os.path.join("work/cna/", name_pattern, "log", name_pattern)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        return {"pipeline_id": "ENSEMBL"}

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8192M",
        )


class cbioportalCnaFilesStepPart(cbioportalExportStepPart):
    """Generate cbioportal continuous (log2), discrete (pseudo-gistic) & segment tables"""

    #: Step name
    name = "cbioportal_cna"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("log2", "gistic")

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = (
            self.config.mapping_tool
            + "."
            + self.config.copy_number_alteration.copy_number_tool
            + ".{library_name}"
        )
        self.input_tpl = os.path.join("work/cna", name_pattern, "out", name_pattern + ".cna")

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        if action == "log2":
            return {
                "action_type": "log2",
                "mappings": self.config.path_gene_id_mappings,
                "extra_args": {"pipeline_id": "ENSEMBL"},
            }
        if action == "gistic":
            return {
                "action_type": "gistic",
                "mappings": self.config.path_gene_id_mappings,
                "extra_args": {
                    "pipeline_id": "ENSEMBL",
                    "amplification": "9",
                },  # Amplification: cn >= 9 (https://doi.org/10.1038/s41586-022-04738-6)
            }

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        return "work/upload/data_cna_{action}.txt".format(action=action)

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8192M",
        )


class cbioportalSegmentStepPart(cbioportalExportStepPart):
    #: Step name
    name = "cbioportal_segment"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_segment.txt"

    def __init__(self, parent):
        super().__init__(parent)
        name_pattern = (
            self.config.mapping_tool
            + "."
            + self.config.copy_number_alteration.copy_number_tool
            + ".{library_name}"
        )
        self.input_tpl = os.path.join(
            self.config.copy_number_alteration.path_copy_number,
            "output",
            name_pattern,
            "out",
            name_pattern + "_dnacopy.seg",
        )

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8192M",
        )


class cbioportalExpressionStepPart(cbioportalExportStepPart):
    """Generates the expression table (z-scores or FPKM)"""

    #: Step name
    name = "cbioportal_expression"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "RNA"

    #: Action
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_expression.txt"

    def __init__(self, parent):
        super().__init__(parent)

        name_pattern = self.config.expression.expression_tool + ".{library_name}"
        self.input_tpl = os.path.join(
            self.config.expression.path_ngs_mapping,
            "output",
            name_pattern,
            "out",
            name_pattern + ".GeneCounts.tab",
        )

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        return {
            "action_type": "expression",
            "extra_args": {
                "pipeline_id": "ENSEMBL",
                "tx_obj": self.parent.w_config.static_data_config.features.path,
            },
        }

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="02:00:00",  # 2 hours
            memory="8192M",
        )


class cbioportalMetaFilesStepPart(BaseStepPart):
    """Generate cbioportal meta data files"""

    #: Step name
    name = "cbioportal_meta_files"

    #: Actions
    actions = ("run",)

    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield from [os.path.join("work/upload", f) for f in META_FILES["always_present"]]
        yield from [os.path.join("work/upload", f) for f in META_FILES["sequenced"]]
        if self.config.copy_number_alteration.enabled:
            yield from [os.path.join("work/upload", f) for f in META_FILES["cna"]]
            if self.config.study.reference_genome == "hg19":
                yield from [os.path.join("work/upload", f) for f in META_FILES["segment"]]
        if self.config.expression.enabled:
            yield from [os.path.join("work/upload", f) for f in META_FILES["rna_seq_mrna"]]

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        return self.config.study.model_dump(by_alias=True)


class cbioportalClinicalDataStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_clinical_data"

    #: Actions
    actions = ("run",)

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        donors = {}
        for extraction_type in ("DNA", "RNA"):
            if extraction_type == "RNA" and not self.config.expression.enabled:
                continue
            self.extraction_type = extraction_type
            for lib in self._yield_libraries():
                assert lib.test_sample.extra_infos["extractionType"] == extraction_type
                donor_name = lib.test_sample.bio_sample.bio_entity.name
                sample_name = lib.test_sample.bio_sample.name
                if donor_name not in donors:
                    donors[donor_name] = {}
                if sample_name not in donors[donor_name]:
                    donors[donor_name][sample_name] = {}
                # Multiple libraries should not be returned by _yield_libraries
                assert extraction_type not in donors[donor_name][sample_name]
                donors[donor_name][sample_name][extraction_type] = lib.name
        assert "__config" not in donors.keys(), "__config is a reserved key, not a valid donor"
        donors["__config"] = self.config.model_dump(by_alias=True)
        # donors["__config"]["vcf2maf"] = dict(self.config.vcf2maf)
        # donors["__config"]["study"] = dict(self.config.study)
        return donors

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        for data_type in ("patient", "sample"):
            yield data_type, "work/upload/data_clinical_{data_type}.txt".format(data_type=data_type)


class cbioportalCaseListsStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_case_lists"

    #: Actions
    actions = ("run",)

    def get_args(self, action):
        # Validate action
        self._validate_action(action)
        samples = dict(zip(CASE_LIST_FILES.keys(), [[] for _ in range(len(CASE_LIST_FILES))]))
        self.extraction_type = "DNA"
        for lib in self._yield_libraries():
            sample_name = lib.test_sample.bio_sample.name
            samples["sequenced"] += [sample_name]
            if self.config.copy_number_alteration.enabled:
                samples["cna"] += [sample_name]
                samples["cnaseq"] += [sample_name]
        if self.config.expression.enabled:
            self.extraction_type = "RNA"
            for lib in self._yield_libraries():
                sample_name = lib.test_sample.bio_sample.name
                samples["rna_seq_mrna"] += [sample_name]
        args = {}
        for k, v in samples.items():
            if len(v) > 0:
                args[k] = CASE_LIST_FILES[k]
                args[k]["samples"] = v
        if self.config.copy_number_alteration.enabled and self.config.expression.enabled:
            args["3way_complete"] = CASE_LIST_FILES["3way_complete"]
            args["3way_complete"]["samples"] = []
            for sample_name in args["cnaseq"]["samples"]:
                if sample_name in args["rna_seq_mrna"]["samples"]:
                    args["3way_complete"]["samples"] += [sample_name]
        args["__cancer_study_id"] = self.config.study.cancer_study_id
        return args

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        case_lists = {}
        case_lists["sequenced"] = CASE_LIST_FILES["sequenced"]["filename"]
        if self.config.copy_number_alteration.enabled:
            case_lists["cna"] = CASE_LIST_FILES["cna"]["filename"]
            case_lists["cnaseq"] = CASE_LIST_FILES["cnaseq"]["filename"]
        if self.config.expression.enabled:
            case_lists["rna_seq_mrna"] = CASE_LIST_FILES["rna_seq_mrna"]["filename"]
            if self.config.copy_number_alteration.enabled:
                case_lists["3way_complete"] = CASE_LIST_FILES["3way_complete"]["filename"]
        for case, filename in case_lists.items():
            yield case, os.path.join("work/upload/case_lists", filename)


# ================================================================================================
#
# cBioPortal step class
#
# ================================================================================================


class cbioportalExportWorkflow(BaseStep):
    """Perform cbioportal preparation"""

    #: Workflow name
    name = "cbioportal_export"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
            config_model_class=CbioportalExportConfigModel,
        )

        # cBioPortal requires the genome release as GRC[hm]3[78] in the MAF file
        # and hg19 or hg38 in the meta study & meta segment files.
        # Note that segmentation visualisation is not yet implemented for hg38 (2023-06-21)
        # For the user's convenience, the configuration is augmented automatically,
        # before the sub steps are registered, so they all have the updated config.
        translated = "unknown"
        if self.config.vcf2maf.ncbi_build == "GRCh37":
            translated = "hg19"
        if self.config.vcf2maf.ncbi_build == "GRCh38":
            translated = "hg38"
        if self.config.vcf2maf.ncbi_build in ("mm9", "mm10", "GRCm37", "GRCm38", "GRCm39"):
            translated = "mouse"
        self.config.study.reference_genome = translated

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                cbioportalMetaFilesStepPart,
                cbioportalCaseListsStepPart,
                cbioportalClinicalDataStepPart,
                cbioportalCns2CnaStepPart,
                cbioportalCnaFilesStepPart,
                cbioportalSegmentStepPart,
                cbioportalMutationsStepPart,
                cbioportalVcf2MafStepPart,
                cbioportalExpressionStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow(
            self.config.somatic_variant_step,
            workdir=self.config.path_somatic_variant,
            sub_workflow_name="somatic_variant",
        )
        if self.config.copy_number_alteration.enabled:
            if self.config.copy_number_alteration.copy_number_tool in (
                "cnvkit",
                "purecn",
                "sequenza",
            ):
                self.register_sub_workflow(
                    "somatic_targeted_seq_cnv_calling",
                    workdir=self.config.copy_number_alteration.path_copy_number,
                    sub_workflow_name="copy_number_step",
                )
            else:
                self.register_sub_workflow(
                    "somatic_wgs_cnv_calling",
                    workdir=self.config.copy_number_alteration.path_copy_number,
                    sub_workflow_name="copy_number_step",
                )
        if self.config.expression.enabled:
            self.register_sub_workflow(
                "ngs_mapping",
                workdir=self.config.expression.path_ngs_mapping,
                sub_workflow_name="ngs_mapping",
            )

    @listify
    def get_result_files(self):
        yield from self.sub_steps["cbioportal_meta_files"].get_output_files("run")
        yield from self.sub_steps["cbioportal_clinical_data"].get_output_files("run").values()
        yield from self.sub_steps["cbioportal_case_lists"].get_output_files("run").values()

        result_files = []
        result_files += [self.sub_steps["cbioportal_mutations"].get_output_files("run")]
        if self.config.copy_number_alteration.enabled:
            result_files += [self.sub_steps["cbioportal_cna"].get_output_files("log2")]
            result_files += [self.sub_steps["cbioportal_cna"].get_output_files("gistic")]
            if self.config.study.reference_genome == "hg19":
                result_files += [self.sub_steps["cbioportal_segment"].get_output_files("run")]
        if self.config.expression.enabled:
            result_files += [self.sub_steps["cbioportal_expression"].get_output_files("run")]

        yield from result_files

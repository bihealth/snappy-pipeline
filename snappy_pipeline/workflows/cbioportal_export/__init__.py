# -*- coding: utf-8 -*-
"""Implementation of the ``cbioportal_export`` step

This step takes as input a wide range of files previously produced by the
snappy_pipeline. It does the necessary data transformations to comply with the
formats expected by cBioPortal and writes out the necessary metadata and data
files.

The output is a archive (tarball) that can be uploaded to a cbioportal instance
for import.
"""

from collections import OrderedDict
import csv
import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

# cbioportal meta data files
META_FILES = (
    "meta_clinical_patient.txt",
    "meta_clinical_sample.txt",
    "meta_CNA_gistic.txt",
    "meta_CNA_log2.txt",
    "meta_expression_zscores.txt",
    "meta_mutation_extended.txt",
    "meta_segment.txt",
)

# files holding clinical data
CLINICAL_DATA_FILES = ("data_clinical_patients.txt", "data_clinical_samples.txt")

# files with CNA results
CNA_DATA_FILES = ("data_CNA_gistic.txt", "data_CNA_log2.txt", "data_segment.txt")


DEFAULT_CONFIG = r"""
step_config:
  cbioportal_export:
    # Paths to snappy steps containing results to be uploaded
    path_ngs_mapping: ../ngs_mapping                                        # REQUIRED
    path_gene_expression_quantification: ../gene_expression_quantification  # Set to '' in absence of mRNA_seq
    path_somatic_variant_filtration: ../somatic_variant_filtration          # REQUIRED
    path_copy_number_step: ../somatic_targeted_seq_cnv_calling              # Set to '' for panels, and change for WGS
    # Select tools & filter set
    cnv_tool: copywriter                                                    # Other option: cnvkit, Control_FREEC (unsupported)
    tools_somatic_variant_calling: [ "mutect" ]                             # Possibly scalpel, mutect2, strelka2
    filter_set: dkfz_only                                                   # Possibly dkfz_and_ebfilter, dkfz_and_ebfilter_and_oxog, ...
    exclude_variant_with_flag: LowFisherScore  # REQUIRED
    # Additional parameters
    vep_data_path: REQUIRED   # Variant Effect Predictor DB for vcf -> maf conversion, must match genome release
    filter_vcf: REQUIRED      # Common germline variants for vcf -> maf conversion, must match genome release
    ncbi_build: GRCh37
    # Description of dataset in cBioPortal
    type_of_cancer: REQUIRED # REQUIRED
    cancer_study_id: REQUIRED # REQUIRED
    study_description: REQUIRED # REQUIRED
    study_name: REQUIRED # REQUIRED
    study_name_short: REQUIRED # REQUIRED
    cache_version: 100

"""


class cbioportalExportStepPart(BaseStepPart):
    """Base class for gene expression quantifiers"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def iterate_over_biomedsheets(self, action):
        """Iterate over biomedsheet, yield input files for cbioportal CNV aggregation
        or for z-score computation (only return samples with matched RNA-seq"""
        copy_number_step = self.parent.sub_workflows["copy_number_step"]
        cnv_tool = self.config["cnv_tool"]
        if action == "get_zscores_input":
            gene_expression = self.parent.sub_workflows.get("gene_expression_quantification")

        cnv_tpl = "output/{mapper}.{tool}.{library_name}/out/{mapper}.{tool}.{library_name}_"
        exp_tpl = "output/{mapper}.{tool}.{library_name}/out/{mapper}.{tool}.{library_name}.tsv"

        for sheet in filter(is_not_background, self.parent.sheets):
            for donor in sheet.bio_entities.values():
                for biosample in donor.bio_samples.values():
                    if biosample.extra_infos["isTumor"]:
                        sample_name = biosample.name
                        rna_library = None
                        dna_library = None
                        for test_sample in biosample.test_samples.values():
                            for lib in test_sample.ngs_libraries.values():
                                # pick one rna and one dna library -
                                # even if there is more than one lib, cbioportal cannot use it
                                if lib.extra_infos["libraryType"] == "mRNA_seq":
                                    rna_library = lib.name
                                elif lib.extra_infos["libraryType"] == "WES":
                                    dna_library = lib.name
                                elif lib.extra_infos["libraryType"] == "WGS":
                                    dna_library = lib.name

                        # for each tumor sample, format file for DNA and optionally RNA
                        cnv_file = cnv_tpl.format(
                            tool=cnv_tool, mapper="bwa", library_name=dna_library
                        )
                        exp_file = None
                        if rna_library:
                            exp_file = exp_tpl.format(
                                tool="featurecounts", mapper="star", library_name=rna_library
                            )

                        # CNVs
                        if action == "gistic":
                            yield copy_number_step(cnv_file + "gene_call.txt")

                        elif action == "log2":
                            yield copy_number_step(cnv_file + "gene_log2.txt")

                        elif action == "segments":
                            yield copy_number_step(cnv_file + "segments.txt")

                        # Z scores
                        elif action == "get_zscores_input" and dna_library and rna_library:
                            cnv_file_ = copy_number_step(cnv_file + "gene_call.txt")
                            exp_file_ = gene_expression(exp_file)
                            yield sample_name, cnv_file_, exp_file_


class CbioportalStudyMetaFilesStepPart(cbioportalExportStepPart):
    """Generate cbioportal study meta data files"""

    #: Step name
    name = "cbioportal_study_meta_files"


class cbioportalMetaFilesStepPart(cbioportalExportStepPart):
    """Generate cbioportal meta data files"""

    #: Step name
    name = "cbioportal_meta_files"

    @listify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield from [os.path.join("work/upload", f) for f in META_FILES]


class cbioportalClinicalDataStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_clinical_data"

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "patients_tsv", os.path.join("work/upload", "data_clinical_patients.txt")
        yield "samples_tsv", os.path.join("work/upload", "data_clinical_samples.txt")

    def get_sample_sheets(self, action):
        sheets = self.parent.sheets
        return [sheet for sheet in sheets]


class cbioportalCaseListsStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_case_lists"

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        yield "sequenced", "work/upload/case_lists/all_cases_with_mutation_data.txt"


class cbioportalVcf2MafStepPart(cbioportalExportStepPart):
    """Helper class for VCF2MAF step"""

    #: Step name
    name = "cbioportal_vcf2maf"

    def get_args(self, action):
        # Validate action
        self._validate_action(action)

        def args_function(wildcards):
            result = {
                "tumor_sample": wildcards.tumor_library,
                "normal_sample": self.get_normal_lib_name(wildcards),
                "tumor_id": self.get_tumor_bio_sample(wildcards),
                "normal_id": self.get_normal_bio_sample(wildcards),
            }
            return result

        return args_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_tumor_bio_sample(self, wildcards):
        """Return bio sample to tumor dna ngs library"""
        return "-".join(wildcards.tumor_library.split("-")[0:2])

    def get_normal_bio_sample(self, wildcards):
        """Return normal bio sample to tumor dna ngs library"""
        normal_lib = self.get_normal_lib_name(wildcards)
        return "-".join(normal_lib.split("-")[0:2])

    @dictify
    def get_output_files(self, action):
        """Return maf output file"""
        # Validate action
        self._validate_action(action)
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}."
            "{filter_set}.{exon_list}"
        )
        yield "maf", os.path.join("work/maf", name_pattern, "out", name_pattern + ".maf")

    @dictify
    def get_input_files(self, action):
        """Return input vcf for each output maf"""
        # Validate action
        self._validate_action(action)
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}."
            "{filter_set}.{exon_list}"
        )
        tpl = (os.path.join("output", name_pattern, "out", name_pattern + ".vcf.gz"),)
        somatic_variant_filtration = self.parent.sub_workflows["somatic_variant_filtration"]
        yield "vcf", somatic_variant_filtration(tpl)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="12:00:00",  # 12 hours
            memory="5120M",
        )


class cbioportalMafStepPart(cbioportalExportStepPart):
    """Helper class to get all required maf files for cbioportal."""

    #: Step name
    name = "cbioportal_maf"

    @listify
    def get_input_files(self, action):
        """Return list of all input files"""
        # Validate action
        self._validate_action(action)
        name_pattern = (
            "{mapper}.{caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}."
            "{filter_set}.{exon_list}"
        )
        filter_sets = self.config["filter_set"]
        exon_lists = ["genome_wide"]
        callers = set(self.config["tools_somatic_variant_calling"])
        yield from self._yield_result_files_matched(
            os.path.join("work/maf", name_pattern, "out", name_pattern + "{ext}"),
            mapper="bwa",
            caller=callers,
            filter_set=filter_sets,
            exon_list=exon_lists,
            ext=".maf",
        )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl, tumor_library=[sample_pair.tumor_sample.dna_ngs_library.name], **kwargs
                )


class cbioportalCnaFilesStepPart(cbioportalExportStepPart):
    """Generate cbioportal cna data files"""

    #: Step name
    name = "cbioportal_cna_data"

    #: Class available actions
    actions = ("log2", "gistic", "segments")

    @listify
    def get_input_files(self, action):
        # Validate action
        self._validate_action(action)
        yield from self.iterate_over_biomedsheets(action)

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="01:00:00",  # 1 hour
            memory="8G",
        )


class cbioportalZscoresStepPart(cbioportalExportStepPart):
    """Generate a dataframe holding each biosample and the location of the CNV
    calling results and the expression quantification files"""

    #:Step name
    name = "cbioportal_zscores"

    def get_input_files(self, action):
        assert action == "get_zscores_input"
        yield from self.iterate_over_biomedsheets(action)

    def get_df(self, output):
        with open(output[0], "w") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(["ID", "gene_call_filename", "count_filename"])
            for sample_name, cnv_file, exp_file in self.iterate_over_biomedsheets(
                "get_zscores_input"
            ):
                writer.writerow([sample_name, cnv_file, exp_file])


class CbioportalComputeZscoresStepPart(cbioportalExportStepPart):
    """Compute Z-Scores"""

    #:Step name
    name = "cbioportal_compute_zscores"


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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                CbioportalStudyMetaFilesStepPart,
                cbioportalMetaFilesStepPart,
                cbioportalCaseListsStepPart,
                cbioportalClinicalDataStepPart,
                cbioportalCnaFilesStepPart,
                cbioportalMafStepPart,
                cbioportalVcf2MafStepPart,
                cbioportalZscoresStepPart,
                CbioportalComputeZscoresStepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        if self.config["path_somatic_variant_filtration"]:
            self.register_sub_workflow(
                "somatic_variant_filtration", self.config["path_somatic_variant_filtration"]
            )
        if self.config["path_gene_expression_quantification"]:
            self.register_sub_workflow(
                "gene_expression_quantification", self.config["path_gene_expression_quantification"]
            )
        if self.config["path_copy_number_step"]:
            if self.config["cnv_tool"] in ["cnvetti_on_target_postprocess", "copywriter"]:
                self.register_sub_workflow(
                    "somatic_targeted_seq_cnv_calling",
                    workdir=self.config["path_copy_number_step"],
                    sub_workflow_name="copy_number_step",
                )
            else:
                self.register_sub_workflow(
                    "somatic_wgs_cnv_calling",
                    workdir=self.config["path_copy_number_step"],
                    sub_workflow_name="copy_number_step",
                )

        if not self.config.get("tools_somatic_variant_calling", ""):
            self.config["tools_somatic_variant_calling"] = self.w_config["step_config"][
                "somatic_variant_calling"
            ]["tools"]

    @listify
    def get_result_files(self):
        result_files = ("meta_study.txt",) + META_FILES + CLINICAL_DATA_FILES
        if self.config["path_somatic_variant_filtration"]:
            result_files += ("data_mutation_extended.txt",)
            result_files += ("case_lists/all_cases_with_mutation_data.txt",)
        if self.config["path_copy_number_step"]:
            result_files += CNA_DATA_FILES
        if self.config["path_gene_expression_quantification"]:
            result_files += ("data_expression_zscores.txt",)

        yield from [os.path.join("work/upload", f) for f in result_files]

    def check_config(self):
        """Check config attributes for presence"""
        if self.config["cnv_tool"] not in [
            "cnvetti_on_target_postprocess",
            "copywriter",
            "control_freec",
        ]:
            raise MissingConfiguration(msg="Please select a supported tool for the CNV calls")

# -*- coding: utf-8 -*-
"""Tests for the cbioportal_export workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.cbioportal_export import cbioportalExportWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          features:
            path: /path/to/features.gtf

        step_config:
          ngs_mapping:
            star:
          cbioportal_export:
            # Paths to snappy steps containing results to be uploaded
            path_ngs_mapping: /NGS_MAPPING
            expression_tool: star
            path_somatic_variant: /SOM_VAR_FILTRATION
            somatic_variant_calling_tool: mutect2
            somatic_variant_annotation_tool: "vep"
            filter_set: dkfz_only
            path_copy_number: /COPY_NUMBER
            copy_number_tool: cnvkit
            exclude_variant_with_flag: LowFisherScore
            # Description of dataset in cBioPortal
            study:
              type_of_cancer: mixed
              cancer_study_id: mixed_pedion_a02p
              study_description: "PeDiOn project A02P"
              study_name: PeDiOn_A02P
              study_name_short: A02P

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: matched_cancer
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def cbioportal_export_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    mocker,
):
    """Return cbioportalExportWorkflow object pre-configured with cancer sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "/NGS_MAPPING/" + x,
        "somatic_variant": lambda x: "/SOM_VAR_FILTRATION/" + x,
        "copy_number_step": lambda x: "/COPY_NUMBER/" + x,
    }
    # Construct the workflow object
    return cbioportalExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for cbioportalMetaFilesStepPart   -----------------------------------------------------


def test_cbioportal_meta_files_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_input_files("cbioportal_meta_files", "run")


def test_cbioportal_meta_files_step_part_get_output_files(cbioportal_export_workflow):
    """Tests CbioportalStudyMetaFilesStepPart.get_log_file()"""
    # Define expected: all meta files as somatic variants, CNA, segmentation & expression are present
    expected = [
        "work/upload/meta_{}.txt".format(x)
        for x in (
            "study",
            "clinical_patient",
            "clinical_sample",
            "mutation_extended",
            "cna_gistic",
            "cna_log2",
            "segment",
            "expression",
        )
    ]
    actual = list(cbioportal_export_workflow.get_output_files("cbioportal_meta_files", "run"))
    assert actual == expected


def test_cbioportal_meta_files_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_meta_files", "run")


def test_cbioportal_meta_files_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_meta_files", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalClinicalDataStepPart ---------------------------------------------------------


def test_cbioportal_clinical_data_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(AssertionError):
        cbioportal_export_workflow.get_input_files("cbioportal_clinical_data", "run")


def test_cbioportal_clinical_data_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_output_files()"""
    # Define expected
    expected = {
        "patient": "work/upload/data_clinical_patient.txt",
        "sample": "work/upload/data_clinical_sample.txt",
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_clinical_data", "run")
    assert actual == expected


def test_cbioportal_clinical_data_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_clinical_data"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_clinical_data", "run")
    assert actual == expected


def test_cbioportal_clinical_data_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource(
            "cbioportal_clinical_data", "run", resource
        )
        assert actual == expected, msg_error


def test_cbioportal_clinical_data_step_part_get_args(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_args()"""
    # Define expected: all patients & samples
    expected = {
        "P001": {"P001-T1": {"DNA": "P001-T1-DNA1-WGS1", "RNA": "P001-T1-RNA1-mRNA_seq1"}},
        "P002": {
            "P002-T1": {"DNA": "P002-T1-DNA1-WGS1"},
            "P002-T2": {"DNA": "P002-T2-DNA1-WGS1", "RNA": "P002-T2-RNA1-mRNA_seq1"},
        },
    }
    actual = cbioportal_export_workflow.get_args("cbioportal_clinical_data", "run")
    assert actual == expected


# Tests for cbioportalCaseListsStepPart   ----------------------------------------------------------


def test_cbioportal_case_lists_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(AssertionError):
        cbioportal_export_workflow.get_input_files("cbioportal_case_lists", "run")


def test_cbioportal_case_lists_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_output_files()"""
    # Define expected
    expected = {
        "sequenced": "work/upload/case_lists/all_cases_with_mutation_data.txt",
        "cna": "work/upload/case_lists/all_cases_with_cna_data.txt",
        "rna_seq_mrna": "work/upload/case_lists/all_cases_with_mrna_rnaseq_data.txt",
        "cnaseq": "work/upload/case_lists/all_cases_with_mutation_and_cna_data.txt",
        "3way_complete": "work/upload/case_lists/all_cases_with_mutation_and_cna_and_mrna_data.txt",
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_case_lists", "run")
    assert actual == expected


def test_cbioportal_case_lists_data_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_case_lists"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_case_lists", "run")
    assert actual == expected


def test_cbioportal_case_lists_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_case_lists", "run", resource)
        assert actual == expected, msg_error


def test_cbioportal_case_lists_step_part_get_args(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_args()"""
    # Define expected: all patients & samples
    expected = {
        "sequenced": {
            "filename": "all_cases_with_mutation_data.txt",
            "name": "Sequenced tumors",
            "description": "Tumors with somatic variant calls",
            "stable_id": "sequenced",
            "category": "all_cases_with_mutation_data",
            "samples": ["P001-T1", "P002-T1", "P002-T2"],
        },
        "cna": {
            "filename": "all_cases_with_cna_data.txt",
            "name": "Tumors with CNA data",
            "description": "Tumors with somatic Copy Number Alteration calls",
            "stable_id": "cna",
            "category": "all_cases_with_cna_data",
            "samples": ["P001-T1", "P002-T1", "P002-T2"],
        },
        "rna_seq_mrna": {
            "filename": "all_cases_with_mrna_rnaseq_data.txt",
            "name": "Tumors with expression data",
            "description": "Tumors with mRNA seq expression data",
            "stable_id": "rna_seq_mrna",
            "category": "all_cases_with_mrna_rnaseq_data",
            "samples": ["P001-T1", "P002-T2"],
        },
        "cnaseq": {
            "filename": "all_cases_with_mutation_and_cna_data.txt",
            "name": "Sequenced tumors with CNA",
            "description": "Tumors with somatic variant & CNA calls",
            "stable_id": "cnaseq",
            "category": "all_cases_with_mutation_and_cna_data",
            "samples": ["P001-T1", "P002-T1", "P002-T2"],
        },
        "3way_complete": {
            "filename": "all_cases_with_mutation_and_cna_and_mrna_data.txt",
            "name": "Sequenced tumors with CNA & expression",
            "description": "Tumors with somatic variant calls, CNA calls & expression data",
            "stable_id": "3way_complete",
            "category": "all_cases_with_mutation_and_cna_and_mrna_data",
            "samples": ["P001-T1", "P002-T2"],
        },
    }
    actual = cbioportal_export_workflow.get_args("cbioportal_case_lists", "run")
    assert actual == expected


# Tests for cbioportalVcf2MafStepPart   ----------------------------------------------------------


def test_cbioportal_vcf2maf_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_input_files()"""
    expected = {
        "vcf": (
            "/SOM_VAR_FILTRATION/output/{mapper}.{caller}.{annotator}."
            "dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}.{exon_list}/out/"
            "{mapper}.{caller}.{annotator}.dkfz_bias_filter.eb_filter."
            "{tumor_library}.{filter_set}.{exon_list}.vcf.gz"
        )
    }
    actual = cbioportal_export_workflow.get_input_files("cbioportal_vcf2maf", "run")
    assert actual == expected


def test_cbioportal_vcf2maf_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_output_files()"""
    # Define expected
    expected = {
        "maf": (
            "work/maf/{mapper}.{caller}.{annotator}.dkfz_bias_filter.eb_filter."
            "{tumor_library}.{filter_set}.{exon_list}/out/{mapper}.{caller}."
            "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
            "{exon_list}.maf"
        )
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_vcf2maf", "run")
    assert actual == expected


def test_cbioportal_vcf2maf_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_log_file()"""
    # Define expected
    base_name_out = (
        "work/maf/{mapper}.{caller}.{annotator}.dkfz_bias_filter.eb_filter."
        "{tumor_library}.{filter_set}.{exon_list}/log/{mapper}.{caller}."
        "{annotator}.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
        "{exon_list}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_vcf2maf", "run")
    assert actual == expected


def test_cbioportal_vcf2maf_step_part_get_args(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_args()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa", "tumor_library": "P001-T1-DNA1-WGS1"})
    expected = {
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "tumor_id": "P001-T1",
        "normal_id": "P001-N1",
    }
    actual = cbioportal_export_workflow.get_args("cbioportal_vcf2maf", "run")(wildcards)
    assert actual == expected


def test_cbioportal_vcf2maf_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 2, "time": "02:00:00", "memory": "5120M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_vcf2maf", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalMutationsStepPart   ----------------------------------------------------------------


def test_cbioportal_mutations_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalMutationsStepPart.get_input_files()"""
    sample = "P00{i}-T{t}"
    base_name = (
        "work/maf/bwa.mutect2.vep.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.dkfz_only.genome_wide/out/"
        "bwa.mutect2.vep.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.dkfz_only.genome_wide.maf"
    )
    expected = {
        sample.format(i=i, t=t): base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))
    }
    actual = cbioportal_export_workflow.get_input_files("cbioportal_mutations", "run")
    assert actual == expected


def test_cbioportal_mutations_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalMutationsStepPart.get_output_files()"""
    # Define expected: default defined workflow.abstract
    expected = "work/upload/data_mutation_extended.txt"
    # Actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_mutations", "run")
    assert actual == expected


def test_cbioportal_mutations_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalMutationsStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_mutations"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_mutations", "run")
    assert actual == expected


def test_cbioportal_mutations_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalMutationsStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_mutations", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalCns2CnaStepPart   -----------------------------------------------------------


def test_cbioportal_cns2cna_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalCns2CnaStepPart.get_input_files()"""
    expected = {
        "DNAcopy": "/COPY_NUMBER/output/{mapper}.{caller}.{tumor_library}/out/{mapper}.{caller}.{tumor_library}_dnacopy.seg"
    }
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cns2cna", "run")
    assert actual == expected


def test_cbioportal_cns2cna_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalCns2CnaStepPart.get_output_files()"""
    # Define expected
    expected = {
        "cna": "work/cna/{mapper}.{caller}.{tumor_library}/out/{mapper}.{caller}.{tumor_library}.cna"
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_cns2cna", "run")
    assert actual == expected


def test_cbioportal_cns2cna_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalCns2CnaStepPart.get_log_file()"""
    # Define expected
    base_name_out = (
        "work/cna/{mapper}.{caller}.{tumor_library}/log/{mapper}.{caller}.{tumor_library}"
    )
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_cns2cna", "run")
    assert actual == expected


def test_cbioportal_cns2cna_step_part_get_args(cbioportal_export_workflow):
    """Tests cbioportalCns2CnaStepPart.get_args()"""
    expected = {
        "features": "/path/to/features.gtf",
        "pipeline_id": "ENSEMBL",
    }
    actual = cbioportal_export_workflow.get_args("cbioportal_cns2cna", "run")
    assert actual == expected


def test_cbioportal_cns2cna_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalCns2CnaStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 2, "time": "02:00:00", "memory": "8192M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_cns2cna", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalCnaFilesStepPart   -----------------------------------------------------------


def test_cbioportal_cna_step_part_get_input_files_log2(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_input_files() - action 'log2'"""
    # Define expected
    sample = "P00{i}-T{t}"
    base_name = "work/cna/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/out/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.cna"
    expected = {
        sample.format(i=i, t=t): base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))
    }
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cna", "log2")
    assert actual == expected


def test_cbioportal_cna_step_part_get_input_files_gistic(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_input_files() - action 'gistic'"""
    # Define expected
    sample = "P00{i}-T{t}"
    base_name = "work/cna/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/out/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1.cna"
    expected = {
        sample.format(i=i, t=t): base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))
    }
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cna", "gistic")
    assert actual == expected


def test_cbioportal_cna_step_part_get_output_files_log2(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_output_files() - action 'log2'"""
    # Define expected
    expected = "work/upload/data_cna_log2.txt"
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_cna", "log2")
    assert actual == expected


def test_cbioportal_cna_step_part_get_output_files_gistic(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_output_files() - action 'gistic'"""
    # Define expected
    expected = "work/upload/data_cna_gistic.txt"
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_cna", "gistic")
    assert actual == expected


def test_cbioportal_cna_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_cna"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_cna", "log2")
    assert actual == expected


def test_cbioportal_cna_step_part_get_args_log2(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_args() -action 'log2'"""
    # Define expected
    expected = {"action_type": "log2", "extra_args": {"pipeline_id": "ENSEMBL"}}
    # Get actual
    actual = cbioportal_export_workflow.get_args("cbioportal_cna", "log2")
    assert actual == expected


def test_cbioportal_cna_step_part_get_args_gistic(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_args() -action 'gistic'"""
    # Define expected
    expected = {
        "action_type": "gistic",
        "extra_args": {"amplification": "9", "pipeline_id": "ENSEMBL"},
    }
    # Get actual
    actual = cbioportal_export_workflow.get_args("cbioportal_cna", "gistic")
    assert actual == expected


def test_cbioportal_cna_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "02:00:00", "memory": "8192M", "partition": "medium"}
    # Evaluate
    all_actions = cbioportal_export_workflow.substep_getattr("cbioportal_cna", "actions")
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = cbioportal_export_workflow.get_resource("cbioportal_cna", action, resource)
            assert actual == expected, msg_error


# Tests for cbioportalSegmentStepPart   --------------------------------------------------------


def test_cbioportal_segment_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalSegmentStepPart.get_input_files()"""
    # Define expected
    sample = "P00{i}-T{t}"
    base_name = "/COPY_NUMBER/output/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1/out/bwa.cnvkit.P00{i}-T{t}-DNA1-WGS1_dnacopy.seg"
    expected = {
        sample.format(i=i, t=t): base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))
    }
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_segment", "run")
    assert actual == expected


def test_cbioportal_segment_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalSegmentStepPart.get_output_files()"""
    # Define expected
    expected = "work/upload/data_segment.txt"
    # Actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_segment", "run")
    assert actual == expected


def test_cbioportal_segment_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalSegmentStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_segment"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_segment", "run")
    assert actual == expected


def test_cbioportal_segment_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalSegmentStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 2, "time": "02:00:00", "memory": "8192M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_segment", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalExpressionStepPart   -----------------------------------------------------


def test_cbioportal_expression_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalExpressionStepPart.get_input_files()"""
    # Define expected
    sample = "P00{i}-T{t}"
    base_name = (
        "/NGS_MAPPING/output/star.P00{i}-T{t}-RNA1-mRNA_seq1/out/"
        "star.P00{i}-T{t}-RNA1-mRNA_seq1.GeneCounts.tab"
    )
    expected = {sample.format(i=i, t=t): base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 2))}
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_expression", "run")
    assert actual == expected


def test_cbioportal_expression_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalExpressionStepPart.get_output_files()"""
    # Define expected
    expected = "work/upload/data_expression.txt"
    # Actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_expression", "run")
    assert actual == expected


def test_cbioportal_expression_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalExpressionStepPart.get_log_file()"""
    # Define expected
    base_name_out = "work/log/cbioportal_expression"
    expected = get_expected_log_files_dict(base_out=base_name_out, extended=False)
    # Get actual
    actual = cbioportal_export_workflow.get_log_file("cbioportal_expression", "run")
    assert actual == expected


def test_cbioportal_expression_step_part_get_args(cbioportal_export_workflow):
    """Tests cbioportalExpressionStepPart.get_args()"""
    expected = {
        "action_type": "expression",
        "extra_args": {"tx_obj": "/path/to/features.gtf", "pipeline_id": "ENSEMBL"},
    }
    actual = cbioportal_export_workflow.get_args("cbioportal_expression", "run")
    assert actual == expected


def test_cbioportal_expression_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalExpressionStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 2, "time": "02:00:00", "memory": "8192M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_expression", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalExportWorkflow   -------------------------------------------------------------


def test_cbioportal_export_workflow(cbioportal_export_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = [
        "cbioportal_case_lists",
        "cbioportal_clinical_data",
        "cbioportal_cna",
        "cbioportal_cns2cna",
        "cbioportal_expression",
        "cbioportal_meta_files",
        "cbioportal_mutations",
        "cbioportal_segment",
        "cbioportal_vcf2maf",
    ]
    actual = list(sorted(cbioportal_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = [
        "work/upload/case_lists/all_cases_with_cna_data.txt",
        "work/upload/case_lists/all_cases_with_mrna_rnaseq_data.txt",
        "work/upload/case_lists/all_cases_with_mutation_and_cna_and_mrna_data.txt",
        "work/upload/case_lists/all_cases_with_mutation_and_cna_data.txt",
        "work/upload/case_lists/all_cases_with_mutation_data.txt",
        "work/upload/data_clinical_patient.txt",
        "work/upload/data_clinical_sample.txt",
        "work/upload/data_cna_gistic.txt",
        "work/upload/data_cna_log2.txt",
        "work/upload/data_expression.txt",
        "work/upload/data_mutation_extended.txt",
        "work/upload/data_segment.txt",
        "work/upload/meta_clinical_patient.txt",
        "work/upload/meta_clinical_sample.txt",
        "work/upload/meta_cna_gistic.txt",
        "work/upload/meta_cna_log2.txt",
        "work/upload/meta_expression.txt",
        "work/upload/meta_mutation_extended.txt",
        "work/upload/meta_segment.txt",
        "work/upload/meta_study.txt",
    ]
    actual = list(sorted(cbioportal_export_workflow.get_result_files()))
    assert actual == expected

    assert cbioportal_export_workflow.check_config() == 0

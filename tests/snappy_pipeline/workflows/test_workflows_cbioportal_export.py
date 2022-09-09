# -*- coding: utf-8 -*-
"""Tests for the cbioportal_export workflow module code"""

from pathlib import Path
import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.cbioportal_export import cbioportalExportWorkflow

from .conftest import patch_module_fs


@pytest.fixture
def cbioportal_zscores_get_df_content():
    """Return contents of Z-Score data frame."""
    p1 = "P001-T1"
    p2 = "P002-T2"
    som_base = (
        "SOM_CNV_CALLING/output/bwa.copywriter.{p}-DNA1-WGS1/out/"
        "bwa.copywriter.{p}-DNA1-WGS1_gene_call.txt"
    )
    gexp_base = (
        "GENE_EXP_QUANTIFICATION/output/star.featurecounts.{p}-RNA1-mRNA_seq1/out/"
        "star.featurecounts.{p}-RNA1-mRNA_seq1.tsv"
    )
    p1_som = som_base.format(p=p1)
    p1_gexp = gexp_base.format(p=p1)
    p2_som = som_base.format(p=p2)
    p2_gexp = gexp_base.format(p=p2)
    return textwrap.dedent(
        f"""
        ID\tgene_call_filename\tcount_filename
        {p1}\t{p1_som}\t{p1_gexp}
        {p2}\t{p2_som}\t{p2_gexp}
        """
    ).lstrip()


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

        cbioportal_export:
          # Paths to snappy steps containing results to be uploaded
          path_ngs_mapping: ../NGS_MAPPING
          path_gene_expression_quantification: ../GENE_EXP_QUANTIFICATION
          path_somatic_variant_filtration: ../SOM_VAR_FILTRATION
          path_copy_number_step: ../SOM_CNV_CALLING
          # Select tools & filter set
          cnv_tool: control_freec
          tools_somatic_variant_calling: [ "mutect2" ]
          filter_set: dkfz_and_ebfilter
          exclude_variant_with_flag: LowFisherScore
          # Additional parameters
          vep_data_path: /path/to/VEP/static_data/GRCh37
          species: homo_sapiens
          # ncbi_build: GRCh37
          filter_vcf: /path/to/ExAC_nonTCGA.r0.3.1.sites.vep.vcf.gz
          cache_version: 100
          vep_custom: ""
          # Description of dataset in cBioPortal
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
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "copy_number_step": lambda x: "SOM_CNV_CALLING/" + x,
        "somatic_variant_filtration": lambda x: "SOM_VAR_FILTRATION/" + x[0],
        "gene_expression_quantification": lambda x: "GENE_EXP_QUANTIFICATION/" + x,
    }
    # Construct the workflow object
    return cbioportalExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CbioportalStudyMetaFilesStepPart   -----------------------------------------------------


def test_cbioportal_study_meta_files_step_part_get_input_files(cbioportal_export_workflow):
    """Tests CbioportalStudyMetaFilesStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_study_meta_files", "run")


def test_cbioportal_study_meta_files_step_part_get_output_files(cbioportal_export_workflow):
    """Tests CbioportalStudyMetaFilesStepPart.get_output_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_study_meta_files", "run")


def test_cbioportal_study_meta_files_step_part_get_log_file(cbioportal_export_workflow):
    """Tests CbioportalStudyMetaFilesStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_study_meta_files", "run")


def test_cbioportal_study_meta_files_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests CbioportalStudyMetaFilesStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource(
            "cbioportal_study_meta_files", "run", resource
        )
        assert actual == expected, msg_error


# Tests for cbioportalMetaFilesStepPart ------------------------------------------------------------


def test_cbioportal_meta_files_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_input_files("cbioportal_meta_files", "run")(None)


def test_cbioportal_meta_files_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_output_files()"""
    expected = [
        "work/upload/meta_clinical_patient.txt",
        "work/upload/meta_clinical_sample.txt",
        "work/upload/meta_CNA_gistic.txt",
        "work/upload/meta_CNA_log2.txt",
        "work/upload/meta_expression_zscores.txt",
        "work/upload/meta_mutation_extended.txt",
        "work/upload/meta_segment.txt",
    ]
    actual = cbioportal_export_workflow.get_output_files("cbioportal_meta_files", "run")
    assert actual == expected


def test_cbioportal_meta_files_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalMetaFilesStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_meta_files", "run")


def test_cbioportal_meta_files_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests ExpansionHunterStepPart.get_resource_usage()"""
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
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_input_files("cbioportal_clinical_data", "run")(None)


def test_cbioportal_clinical_data_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_output_files()"""
    # Define expected
    expected = {
        "patients_tsv": "work/upload/data_clinical_patients.txt",
        "samples_tsv": "work/upload/data_clinical_samples.txt",
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_clinical_data", "run")
    assert actual == expected


def test_cbioportal_clinical_data_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalClinicalDataStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_clinical_data", "run")


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


# Tests for cbioportalCaseListsStepPart   ----------------------------------------------------------


def test_cbioportal_case_lists_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_input_files("cbioportal_case_lists", "run")(None)


def test_cbioportal_case_lists_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_output_files()"""
    # Define expected
    expected = {"sequenced": "work/upload/case_lists/all_cases_with_mutation_data.txt"}
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_case_lists", "run")
    assert actual == expected


def test_cbioportal_case_lists_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_case_lists", "run")


def test_cbioportal_case_lists_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalCaseListsStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_case_lists", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalVcf2MafStepPart   ----------------------------------------------------------


def test_cbioportal_vcf2maf_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_input_files()"""
    expected = {
        "vcf": (
            "SOM_VAR_FILTRATION/output/{mapper}.{caller}.jannovar_annotate_somatic_vcf."
            "dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}.{exon_list}/out/"
            "{mapper}.{caller}.jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter."
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
            "work/maf/{mapper}.{caller}.jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter."
            "{tumor_library}.{filter_set}.{exon_list}/out/{mapper}.{caller}."
            "jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter.{tumor_library}.{filter_set}."
            "{exon_list}.maf"
        )
    }
    # Get actual
    actual = cbioportal_export_workflow.get_output_files("cbioportal_vcf2maf", "run")
    assert actual == expected


def test_cbioportal_vcf2maf_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalVcf2MafStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_vcf2maf", "run")


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
    expected_dict = {"threads": 4, "time": "12:00:00", "memory": "5120M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_vcf2maf", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalMafStepPart   ----------------------------------------------------------------


def test_cbioportal_maf_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalMafStepPart.get_input_files()"""
    base_name = (
        "work/maf/bwa.mutect.jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.dkfz_only.genome_wide/out/"
        "bwa.mutect.jannovar_annotate_somatic_vcf.dkfz_bias_filter.eb_filter."
        "P00{i}-T{t}-DNA1-WGS1.dkfz_only.genome_wide.maf"
    )
    expected = [base_name.format(i=i, t=t) for i, t in ((1, 1), (2, 1), (2, 2))]
    actual = cbioportal_export_workflow.get_input_files("cbioportal_maf", "run")
    assert actual == expected


def test_cbioportal_maf_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalMafStepPart.get_output_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_maf", "run")


def test_cbioportal_maf_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalMafStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_maf", "run")


def test_cbioportal_maf_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalMafStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_maf", "run", resource)
        assert actual == expected, msg_error


# Tests for cbioportalCnaFilesStepPart   -----------------------------------------------------------


def test_cbioportal_cna_data_step_part_get_input_files_log2(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_input_files() - action 'log2'"""
    # Define expected
    base_name = (
        "SOM_CNV_CALLING/output/bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}/out/"
        "bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}_gene_log2.txt"
    )
    expected = [base_name.format(i=i, t=t, g=g) for i, t, g in ((1, 1, 1), (2, 1, 2), (2, 2, 1))]
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cna_data", "log2")
    assert actual == expected


def test_cbioportal_cna_data_step_part_get_input_files_gistic(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_input_files() - action 'gistic'"""
    # Define expected
    base_name = (
        "SOM_CNV_CALLING/output/bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}/out/"
        "bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}_gene_call.txt"
    )
    expected = [base_name.format(i=i, t=t, g=g) for i, t, g in ((1, 1, 1), (2, 1, 2), (2, 2, 1))]
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cna_data", "gistic")
    assert actual == expected


def test_cbioportal_cna_data_step_part_get_input_files_segments(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_input_files() - action 'segments'"""
    # Define expected
    base_name = (
        "SOM_CNV_CALLING/output/bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}/out/"
        "bwa.copywriter.P00{i}-T{t}-DNA1-WGS{g}_segments.txt"
    )
    expected = [base_name.format(i=i, t=t, g=g) for i, t, g in ((1, 1, 1), (2, 1, 2), (2, 2, 1))]
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_cna_data", "segments")
    assert actual == expected


def test_cbioportal_cna_data_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_output_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_cna_data", "action")


def test_cbioportal_cna_data_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_cna_data", "action")


def test_cbioportal_cna_data_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalCnaFilesStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "01:00:00", "memory": "8G", "partition": "medium"}
    # Evaluate
    all_actions = cbioportal_export_workflow.substep_getattr("cbioportal_cna_data", "actions")
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = cbioportal_export_workflow.get_resource(
                "cbioportal_cna_data", action, resource
            )
            assert actual == expected, msg_error


# Tests for cbioportalZscoresStepPart   ------------------------------------------------------------


def test_cbioportal_zscores_step_part_get_input_files(cbioportal_export_workflow):
    """Tests cbioportalZscoresStepPart.get_input_files()"""
    # Define expected
    p1 = "P001-T1"
    p2 = "P002-T2"
    som_base = (
        "SOM_CNV_CALLING/output/bwa.copywriter.{p}-DNA1-WGS1/out/"
        "bwa.copywriter.{p}-DNA1-WGS1_gene_call.txt"
    )
    gexp_base = (
        "GENE_EXP_QUANTIFICATION/output/star.featurecounts.{p}-RNA1-mRNA_seq1/out/"
        "star.featurecounts.{p}-RNA1-mRNA_seq1.tsv"
    )
    expected = [
        (p1, som_base.format(p=p1), gexp_base.format(p=p1)),
        (p2, som_base.format(p=p2), gexp_base.format(p=p2)),
    ]
    # Get actual
    actual = cbioportal_export_workflow.get_input_files("cbioportal_zscores", "get_zscores_input")
    actual = sorted(actual)
    assert actual == expected


def test_cbioportal_zscores_step_part_get_output_files(cbioportal_export_workflow):
    """Tests cbioportalZscoresStepPart.get_output_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_zscores", "run")


def test_cbioportal_zscores_step_part_get_log_file(cbioportal_export_workflow):
    """Tests cbioportalZscoresStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_zscores", "run")


def test_cbioportal_zscores_step_part_get_df(
    cbioportal_export_workflow, cbioportal_zscores_get_df_content, tmpdir
):
    """Tests cbioportalZscoresStepPart.get_df()"""
    # Define expected
    expected = cbioportal_zscores_get_df_content
    # Prepare method input in temporary directory
    _tmp_work_dir = tmpdir.mkdir("work")
    output_path = str(_tmp_work_dir.join("zscores_mapping_df.tsv"))
    output = (output_path,)
    # Get actual and assert
    cbioportal_export_workflow.substep_getattr("cbioportal_zscores", "get_df")(output)
    actual = Path(output_path).read_text(encoding="utf8")
    assert actual == expected


def test_cbioportal_zscores_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests cbioportalZscoresStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource("cbioportal_zscores", "run", resource)
        assert actual == expected, msg_error


# Tests for CbioportalComputeZscoresStepPart   -----------------------------------------------------


def test_cbioportal_compute_zscores_step_part_get_input_files(cbioportal_export_workflow):
    """Tests CbioportalComputeZscoresStepPart.get_input_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_compute_zscores", "run")


def test_cbioportal_compute_zscores_step_part_get_output_files(cbioportal_export_workflow):
    """Tests CbioportalComputeZscoresStepPart.get_output_files()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_output_files("cbioportal_compute_zscores", "run")


def test_cbioportal_compute_zscores_step_part_get_log_file(cbioportal_export_workflow):
    """Tests CbioportalComputeZscoresStepPart.get_log_file()"""
    # Method not implemented
    with pytest.raises(NotImplementedError):
        cbioportal_export_workflow.get_log_file("cbioportal_compute_zscores", "run")


def test_cbioportal_compute_zscores_step_part_get_resource_usage(cbioportal_export_workflow):
    """Tests CbioportalComputeZscoresStepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    expected_dict = {"threads": 1, "time": "01:00:00", "memory": "2G", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = cbioportal_export_workflow.get_resource(
            "cbioportal_compute_zscores", "run", resource
        )
        assert actual == expected, msg_error


# Tests for cbioportalExportWorkflow   -------------------------------------------------------------


def test_cbioportal_export_workflow(cbioportal_export_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = [
        "cbioportal_case_lists",
        "cbioportal_clinical_data",
        "cbioportal_cna_data",
        "cbioportal_compute_zscores",
        "cbioportal_maf",
        "cbioportal_meta_files",
        "cbioportal_study_meta_files",
        "cbioportal_vcf2maf",
        "cbioportal_zscores",
        "link_out",
    ]
    actual = list(sorted(cbioportal_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = [
        "work/upload/case_lists/all_cases_with_mutation_data.txt",
        "work/upload/data_CNA_gistic.txt",
        "work/upload/data_CNA_log2.txt",
        "work/upload/data_clinical_patients.txt",
        "work/upload/data_clinical_samples.txt",
        "work/upload/data_expression_zscores.txt",
        "work/upload/data_mutation_extended.txt",
        "work/upload/data_segment.txt",
        "work/upload/meta_CNA_gistic.txt",
        "work/upload/meta_CNA_log2.txt",
        "work/upload/meta_clinical_patient.txt",
        "work/upload/meta_clinical_sample.txt",
        "work/upload/meta_expression_zscores.txt",
        "work/upload/meta_mutation_extended.txt",
        "work/upload/meta_segment.txt",
        "work/upload/meta_study.txt",
    ]
    actual = list(sorted(cbioportal_export_workflow.get_result_files()))
    assert actual == expected

# -*- coding: utf-8 -*-
"""Tests for the variant_annotation workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.varfish_export import VarfishExportWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (germline) configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa
            target_coverage_report:
              path_target_interval_list_mapping:
              - name: "Agilent SureSelect Human All Exon V6"
                pattern: "Agilent SureSelect Human All Exon V6*"
                path: "path/to/targets.bed"

          variant_calling:
            tools:
            - gatk3_hc
            gatk3_hc: {}
          variant_annotation:
            tools:
              - vep
            vep: {}
          varfish_export:
            path_ngs_mapping: ../ngs_mapping
            path_variant_calling: ../variant_calling
            path_exon_bed: /path/to/exons.bed
            path_mehari_db: /path/to/mehari.db

        data_sets:
          first_batch:
            file: sheet.tsv
            search_patterns:
            - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
            search_paths: ['/path']
            type: germline_variants
            naming_scheme: only_secondary_id
        """
        ).lstrip()
    )


@pytest.fixture
def varfish_export_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return VarfishExportWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    # Patch out files for aligner indices
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "variant_calling": lambda x: "VAR_CALLING/" + x,
    }
    # Construct the workflow object
    return VarfishExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for MehariAnnotateStepPart -------------------------------------------------------


def test_mehari_step_part_get_input_files_annotate(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_input_files_annotate()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk3_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    # Define expected
    base_name = (
        "VAR_CALLING/output/bwa.gatk3_hc.P001-N1-DNA1-WGS1/out/bwa.gatk3_hc.P001-N1-DNA1-WGS1"
    )
    expected = {
        "ped": "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped",
        "vcf": [base_name + ".vcf.gz"],
    }
    # Get actual
    actual = varfish_export_workflow.get_input_files("mehari", "annotate_seqvars")(wildcards)
    assert actual == expected


def test_mehari_step_part_get_input_files_bam_qc(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_input_files_bam_qc()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "var_caller": "gatk3_hc",
            "index_ngs_library": "P001-N1-DNA1-WGS1",
        }
    )
    # Define expected
    donor_indices = (1, 2, 3)
    base_name_bam = (
        "NGS_MAPPING/output/bwa.P00{i}-N1-DNA1-WGS1/report/bam_qc/bwa.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    base_name_cov = "NGS_MAPPING/output/bwa.P00{i}-N1-DNA1-WGS1/report/alfred_qc/bwa.P00{i}-N1-DNA1-WGS1.alfred.json.gz"
    expected = {
        "bamstats": [base_name_bam.format(i=i, ext="bam.bamstats.txt") for i in donor_indices],
        "flagstats": [base_name_bam.format(i=i, ext="bam.flagstats.txt") for i in donor_indices],
        "idxstats": [base_name_bam.format(i=i, ext="bam.idxstats.txt") for i in donor_indices],
        "alfred_qc": [base_name_cov.format(i=i) for i in donor_indices],
    }
    # Get actual
    actual = varfish_export_workflow.get_input_files("mehari", "bam_qc")(wildcards)
    assert actual == expected


def test_mehari_step_part_get_output_files_annotate(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_output_files_annotate()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.varfish_export.{index_ngs_library}/out/"
        "{mapper}.mehari_annotate_seqvars.{index_ngs_library}"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
        "ped": base_name_out + ".ped",
        "ped_md5": base_name_out + ".ped.md5",
        "output_links": [
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.ped",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.ped.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.gts.tsv.gz",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.gts.tsv.gz.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.db-infos.tsv.gz",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.db-infos.tsv.gz.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.wrapper.py",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.wrapper.py.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.log",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.log.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.conda_info.txt",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.conda_info.txt.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.conda_list.txt",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.conda_list.txt.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.environment.yaml",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_annotate_seqvars.{index_ngs_library}.environment.yaml.md5",
        ],
    }
    # Get actual
    actual = varfish_export_workflow.get_output_files("mehari", "annotate_seqvars")
    assert actual == expected


def test_mehari_step_part_get_output_files_bam_qc(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_output_files_bam_qc()"""
    # Define expected
    base_name_out = (
        "work/{mapper}.varfish_export.{index_ngs_library}/out/"
        "{mapper}.mehari_bam_qc.{index_ngs_library}"
    )
    expected = {
        "bam_qc": base_name_out + ".bam-qc.tsv.gz",
        "bam_qc_md5": base_name_out + ".bam-qc.tsv.gz.md5",
        "output_links": [
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_bam_qc.{index_ngs_library}.bam-qc.tsv.gz",
            "output/{mapper}.varfish_export.{index_ngs_library}/out/{mapper}.mehari_bam_qc.{index_ngs_library}.bam-qc.tsv.gz.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.wrapper.py",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.wrapper.py.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.log",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.log.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.conda_info.txt",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.conda_info.txt.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.conda_list.txt",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.conda_list.txt.md5",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.environment.yaml",
            "output/{mapper}.varfish_export.{index_ngs_library}/log/{mapper}.mehari_bam_qc.{index_ngs_library}.environment.yaml.md5",
        ],
    }
    # Get actual
    actual = varfish_export_workflow.get_output_files("mehari", "bam_qc")
    assert actual == expected


def test_mehari_step_part_get_log_file_annotate(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_log_file()_annotate"""
    # Define expected
    base_name_log = (
        "work/{mapper}.varfish_export.{index_ngs_library}/log/"
        "{mapper}.mehari_annotate_seqvars.{index_ngs_library}"
    )
    base_name_wrapper = (
        "work/{mapper}.varfish_export.{index_ngs_library}/log/"
        "{mapper}.mehari_annotate_seqvars.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
        "env_yaml": base_name_wrapper + ".environment.yaml",
        "env_yaml_md5": base_name_wrapper + ".environment.yaml.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**wrapper_dict, **log_dict}
    # Get actual
    actual = varfish_export_workflow.get_log_file("mehari", "annotate_seqvars")
    assert actual == expected


def test_mehari_step_part_get_log_file_bam_qc(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_log_file()_bam_qc"""
    # Define expected
    base_name_log = (
        "work/{mapper}.varfish_export.{index_ngs_library}/log/"
        "{mapper}.mehari_bam_qc.{index_ngs_library}"
    )
    base_name_wrapper = (
        "work/{mapper}.varfish_export.{index_ngs_library}/log/"
        "{mapper}.mehari_bam_qc.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name_wrapper + ".wrapper.py",
        "wrapper_md5": base_name_wrapper + ".wrapper.py.md5",
        "env_yaml": base_name_wrapper + ".environment.yaml",
        "env_yaml_md5": base_name_wrapper + ".environment.yaml.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name_log)
    expected = {**wrapper_dict, **log_dict}
    # Get actual
    actual = varfish_export_workflow.get_log_file("mehari", "bam_qc")
    assert actual == expected


def test_mehari_step_part_get_params_annotate(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_params_annotate()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {"step_name": "varfish_export"}
    actual = varfish_export_workflow.get_params("mehari", "annotate_seqvars")(wildcards)
    assert actual == expected


def test_mehari_step_part_get_params_bam_qc(varfish_export_workflow):
    """Tests MehariAnnotateStepPart._get_params_bam_qc()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {
        "P001-N1-DNA1-WGS1": "P001-N1-DNA1-WGS1",
        "P002-N1-DNA1-WGS1": "P002-N1-DNA1-WGS1",
        "P003-N1-DNA1-WGS1": "P003-N1-DNA1-WGS1",
    }
    actual = varfish_export_workflow.get_params("mehari", "bam_qc")(wildcards)
    assert actual == expected


def test_mehari_step_part_get_resource_usage(varfish_export_workflow):
    """Tests MehariAnnotateStepPart.get_resource_usage()"""
    all_actions = varfish_export_workflow.substep_getattr("mehari", "actions")
    # Define expected
    expected_dict = {"threads": 2, "time": "1-00:00:00", "memory": "14G", "partition": "medium"}
    # Evaluate
    for action in all_actions:
        for resource, expected in expected_dict.items():
            msg_error = f"Assertion error for resource '{resource}' in action '{action}'."
            actual = varfish_export_workflow.get_resource("mehari", action, resource)()
            assert actual == expected, msg_error


# Tests for VarfishExportWorkflow  ----------------------------------------------------------------


def test_varfish_export_workflow(varfish_export_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "mehari", "write_pedigree"]
    actual = list(sorted(varfish_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/bwa.varfish_export.P00{i}-N1-DNA1-WGS1/{dir_}/"
        "bwa.mehari_{action}.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    # Files in `out` directory
    expected = [
        tpl.format(dir_="out", i=i, ext=ext, action=action)
        for i in (1, 4)  # only for indices
        for (action, ext) in (
            ("annotate_seqvars", "gts.tsv.gz"),
            ("annotate_seqvars", "gts.tsv.gz.md5"),
            ("annotate_seqvars", "db-infos.tsv.gz"),
            ("annotate_seqvars", "db-infos.tsv.gz.md5"),
            ("bam_qc", "bam-qc.tsv.gz"),
            ("bam_qc", "bam-qc.tsv.gz.md5"),
            ("annotate_seqvars", "ped"),
            ("annotate_seqvars", "ped.md5"),
        )
    ]
    # Files in `log` directory
    expected += [
        tpl.format(dir_="log", i=i, ext=ext, action=action)
        for i in (1, 4)  # only for indices
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
            "wrapper.py",
            "wrapper.py.md5",
            "environment.yaml",
            "environment.yaml.md5",
        )
        for action in ("annotate_seqvars", "bam_qc")
    ]
    expected = sorted(expected)
    actual = sorted(varfish_export_workflow.get_result_files())
    assert actual == expected

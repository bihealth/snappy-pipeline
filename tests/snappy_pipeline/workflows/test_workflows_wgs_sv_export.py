# -*- coding: utf-8 -*-
"""Tests for the wgs_sv_export workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.wgs_sv_export import WgsSvExportWorkflow

from .common import get_expected_log_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for germline configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          dbsnp:
            path: /path/to/dbsnp.vcf.gz

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            compute_coverage_bed: true
            path_target_regions: /path/to/regions.bed
            bwa:
              path_index: /path/to/bwa/index.fa

          wgs_sv_export:
            path_wgs_sv_annotation: ../WGS_SV_ANNOTATION
            path_wgs_sv_export: ../WGS_SV_EXPORT
            tools_ngs_mapping: [bwa]
            tools_wgs_sv_calling: [delly2, popdel]

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
def wgs_sv_export_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    germline_sheet_fake_fs,
    mocker,
):
    """Return WgsSvExportWorkflow object pre-configured with germline sheet"""
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", germline_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.wgs_sv_calling", germline_sheet_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really a NGSMappingPipelineStep here
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "wgs_sv_annotation": lambda x: "WGS_SV_ANNOTATION/" + x,
        "wgs_sv_calling": lambda x: "WGS_SV_CALLING/" + x,
    }
    # Construct the workflow object
    return WgsSvExportWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for VarfishAnnotatorAnnotateStepPart  ------------------------------------------------------


def test_varfish_annotator_step_part_call_get_input_files_delly2(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_input_files() with delly2"""
    wildcards = Wildcards(
        fromdict={
            "index_ngs_library": "INDEX_NGS_LIBRARY",
            "mapper": "bwa",
            "var_caller": "delly2",
        }
    )
    wgs_base_name = (
        "WGS_SV_CALLING/output/bwa.delly2.INDEX_NGS_LIBRARY/out/" "bwa.delly2.INDEX_NGS_LIBRARY"
    )
    expected = {
        "ped": "work/write_pedigree.INDEX_NGS_LIBRARY/out/INDEX_NGS_LIBRARY.ped",
        "vcf": wgs_base_name + ".vcf.gz",
        "vcf_md5": wgs_base_name + ".vcf.gz.md5",
        "tbi": wgs_base_name + ".vcf.gz.tbi",
        "tbi_md5": wgs_base_name + ".vcf.gz.tbi.md5",
    }
    actual = wgs_sv_export_workflow.get_input_files("varfish_annotator", "annotate")(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_output_files()"""
    base_name_out = (
        "work/bwa.delly2.varfish_annotated.INDEX_NGS_LIBRARY/out/"
        "bwa.delly2.varfish_annotated.INDEX_NGS_LIBRARY"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "feature_effects": base_name_out + ".feature-effects.tsv.gz",
        "feature_effects_md5": base_name_out + ".feature-effects.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
    }
    actual = wgs_sv_export_workflow.get_output_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_call_get_input_files_popdel(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_input_files() with popdel"""
    wildcards = Wildcards(
        fromdict={
            "index_ngs_library": "INDEX_NGS_LIBRARY",
            "mapper": "bwa",
            "var_caller": "popdel",
        }
    )
    wgs_base_name = (
        "WGS_SV_ANNOTATION/output/bwa.popdel.annotated.INDEX_NGS_LIBRARY/out/"
        "bwa.popdel.annotated.INDEX_NGS_LIBRARY"
    )
    expected = {
        "ped": "work/write_pedigree.INDEX_NGS_LIBRARY/out/INDEX_NGS_LIBRARY.ped",
        "vcf": wgs_base_name + ".vcf.gz",
        "vcf_md5": wgs_base_name + ".vcf.gz.md5",
        "tbi": wgs_base_name + ".vcf.gz.tbi",
        "tbi_md5": wgs_base_name + ".vcf.gz.tbi.md5",
    }
    actual = wgs_sv_export_workflow.get_input_files("varfish_annotator", "annotate")(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_call_get_output_files(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_output_files()"""
    base_name_out = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/out/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    expected = {
        "gts": base_name_out + ".gts.tsv.gz",
        "gts_md5": base_name_out + ".gts.tsv.gz.md5",
        "feature_effects": base_name_out + ".feature-effects.tsv.gz",
        "feature_effects_md5": base_name_out + ".feature-effects.tsv.gz.md5",
        "db_infos": base_name_out + ".db-infos.tsv.gz",
        "db_infos_md5": base_name_out + ".db-infos.tsv.gz.md5",
    }
    actual = wgs_sv_export_workflow.get_output_files("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_call_get_log_file(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_log_file()"""
    # Define expected
    base_name = (
        "work/{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}/log/"
        "{mapper}.{var_caller}.varfish_annotated.{index_ngs_library}"
    )
    wrapper_dict = {
        "wrapper": base_name + ".wrapper.py",
        "wrapper_md5": base_name + ".wrapper.py.md5",
    }
    log_dict = get_expected_log_files_dict(base_out=base_name)
    expected = {**wrapper_dict, **log_dict}
    # Get actual
    actual = wgs_sv_export_workflow.get_log_file("varfish_annotator", "annotate")
    assert actual == expected


def test_varfish_annotator_step_part_get_params(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_params()"""
    wildcards = Wildcards(fromdict={"index_ngs_library": "P001-N1-DNA1-WGS1"})
    expected = {"is_wgs": True, "step_name": "wgs_sv_export", "varfish_server_compatibility": False}
    actual = wgs_sv_export_workflow.get_params("varfish_annotator", "annotate")(wildcards)
    assert actual == expected


def test_varfish_annotator_step_part_get_resource_usage(wgs_sv_export_workflow):
    """Tests VarfishAnnotatorAnnotateStepPart.get_resource_usage()"""
    # Define expected
    expected_dict = {"threads": 2, "time": "4-04:00:00", "memory": "14336M", "partition": "medium"}
    # Evaluate
    for resource, expected in expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}'."
        actual = wgs_sv_export_workflow.get_resource("varfish_annotator", "annotate", resource)
        assert actual == expected, msg_error


# Tests for WgsSvExportWorkflow   ------------------------------------------------------------------


def test_wgs_sv_annotation_workflow(wgs_sv_export_workflow):
    """Tests simple functionality of the workflow."""
    # Check created sub steps
    expected = ["link_out", "varfish_annotator", "write_pedigree"]
    actual = list(sorted(wgs_sv_export_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = (
        "output/{mapper}.{cnv_caller}.varfish_annotated.P00{i}-N1-DNA1-WGS1/{dir_}/"
        "{mapper}.{cnv_caller}.varfish_annotated.P00{i}-N1-DNA1-WGS1.{ext}"
    )
    # Expected files in `out`
    expected = [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, dir_="out", ext=ext)
        for i in (1, 4)
        for ext in (
            "gts.tsv.gz",
            "gts.tsv.gz.md5",
            "db-infos.tsv.gz",
            "db-infos.tsv.gz.md5",
            "feature-effects.tsv.gz",
            "feature-effects.tsv.gz.md5",
        )
        for mapper in ("bwa",)
        for cnv_caller in ("delly2", "popdel")
    ]
    # Expected files in `log`
    expected += [
        tpl.format(mapper=mapper, cnv_caller=cnv_caller, i=i, dir_="log", ext=ext)
        for i in (1, 4)
        for ext in (
            "log",
            "log.md5",
            "conda_info.txt",
            "conda_info.txt.md5",
            "conda_list.txt",
            "conda_list.txt.md5",
        )
        for mapper in ("bwa",)
        for cnv_caller in ("delly2", "popdel")
    ]
    expected = list(sorted(expected))
    actual = list(sorted(wgs_sv_export_workflow.get_result_files()))
    assert actual == expected

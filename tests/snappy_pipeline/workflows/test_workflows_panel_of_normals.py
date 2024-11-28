# -*- coding: utf-8 -*-
"""Tests for the panel_of_normals workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.panel_of_normals import PanelOfNormalsWorkflow

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


@pytest.fixture(scope="module")  # otherwise: performance issues
def minimal_config():
    """Return YAML parsing result for (cancer) configuration"""
    yaml = ruamel_yaml.YAML()
    return yaml.load(
        textwrap.dedent(
            r"""
        static_data_config:
          reference:
            path: /path/to/ref.fa
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz
          features:
            path: /path/to/annotations.gtf

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
            bwa:
              path_index: /path/to/bwa/index.fa

          panel_of_normals:
              tools: ['mutect2', 'cnvkit', 'purecn']
              ignore_chroms: [GL*]
              path_ngs_mapping: ../ngs_mapping
              mutect2:
                  germline_resource: /path/to/germline_resource.vcf
                  path_normals_list: ""
              cnvkit:
                  ignore_chroms: [MT]
                  path_target_interval_list_mapping: []
                  path_normals_list: ""
                  diploid_parx_genome: GRCh38
              purecn:
                  path_normals_list: ""
                  path_target_interval_list_mapping:
                  - name: panel
                    pattern: panel
                    path: /path/to/baits.bed
                  path_genomicsDB: /path/to/mutect2/genomicsDB
                  genome_name: "unknown"

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
def panel_of_normals_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    autobin_result_pon_fake_fs,
    mocker,
):
    """Return PanelOfNormalsWorkflow object pre-configured with germline sheet"""
    # Patch out file-system to enable reading autobin output
    autobin_result_pon_fake_fs.fs.create_file(
        file_path="work/bwa.cnvkit/out/bwa.cnvkit.autobin.txt",
        contents="Target: -1 2000\n",
        create_missing_dirs=True,
    )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.panel_of_normals", autobin_result_pon_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {"ngs_mapping": lambda x: "NGS_MAPPING/" + x}
    # Construct the workflow object
    return PanelOfNormalsWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for Mutect2StepPart ------------------------------------------------------------------------


def test_mutect2_step_part_get_input_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_input_files_prepare_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "normal_library": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "normal_bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "normal_bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "prepare_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
        }
    )
    expected = {
        "normals": [
            "work/bwa.mutect2/out/bwa.mutect2.P001-N1-DNA1-WGS1.prepare.vcf.gz",
            "work/bwa.mutect2/out/bwa.mutect2.P002-N1-DNA1-WGS1.prepare.vcf.gz",
        ],
    }
    actual = panel_of_normals_workflow.get_input_files("mutect2", "create_panel")(wildcards)
    assert actual == expected


def test_mutect2_step_part_get_output_files_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_prepare_panel()"""
    expected = {
        "vcf": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz",
        "vcf_md5": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.md5",
        "vcf_tbi": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.tbi",
        "vcf_tbi_md5": "work/{mapper}.mutect2/out/{mapper}.mutect2.{normal_library}.prepare.vcf.gz.tbi.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_output_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_output_vcf_files_dict(base_out=base_name_out)
    base_name_out = "work/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB"
    expected["db"] = base_name_out + ".tar.gz"
    expected["db_md5"] = expected["db"] + ".md5"
    actual = panel_of_normals_workflow.get_output_files("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_prepare_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_prepare_panel()"""
    base_name_out = "work/{mapper}.mutect2/log/{mapper}.mutect2.{normal_library}.prepare"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "prepare_panel")
    assert actual == expected


def test_mutect2_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests Mutect2StepPart._get_log_files_create_panel()"""
    base_name_out = "work/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals"
    expected = get_expected_log_files_dict(base_out=base_name_out)
    actual = panel_of_normals_workflow.get_log_file("mutect2", "create_panel")
    assert actual == expected


def test_mutect2_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests Mutect2StepPart.get_resource_usage()"""
    # Define expected: default defined workflow.abstract
    create_panel_expected_dict = {
        "threads": 2,
        "time": "48:00:00",
        "memory": "30G",
        "partition": "medium",
    }
    prepare_panel_expected_dict = {
        "threads": 2,
        "time": "3-00:00:00",
        "memory": "8G",
        "partition": "medium",
    }

    # Evaluate action `create_panel`
    for resource, expected in create_panel_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'create_panel'."
        actual = panel_of_normals_workflow.get_resource("mutect2", "create_panel", resource)()
        assert actual == expected, msg_error

    # Evaluate action `prepare_panel`
    for resource, expected in prepare_panel_expected_dict.items():
        msg_error = f"Assertion error for resource '{resource}' for action 'prepare_panel'."
        actual = panel_of_normals_workflow.get_resource("mutect2", "prepare_panel", resource)()
        assert actual == expected, msg_error


# Tests for PureCnStepPart -------------------------------------------------------------------------


def test_purecn_step_part_get_output_files_install(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_install()"""
    expected = {"container": "work/containers/out/purecn.simg"}
    actual = panel_of_normals_workflow.get_output_files("purecn", "install")
    assert actual == expected


def test_purecn_step_part_get_log_file_install(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_install()"""
    expected = get_expected_log_files_dict(base_out="work/containers/log/purecn")
    actual = panel_of_normals_workflow.get_log_file("purecn", "install")
    assert actual == expected


def test_purecn_step_part_get_input_files_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_prepare()"""
    expected = {"container": "work/containers/out/purecn.simg"}
    actual = panel_of_normals_workflow.get_input_files("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_output_files_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_prepare()"""
    expected = {
        "intervals": "work/purecn/out/panel_unknown.list",
        "optimized": "work/purecn/out/panel_unknown.bed.gz",
        "tbi": "work/purecn/out/panel_unknown.bed.gz.tbi",
        "intervals_md5": "work/purecn/out/panel_unknown.list.md5",
        "optimized_md5": "work/purecn/out/panel_unknown.bed.gz.md5",
        "tbi_md5": "work/purecn/out/panel_unknown.bed.gz.tbi.md5",
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_log_file_prepare(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_prepare()"""
    expected = get_expected_log_files_dict(base_out="work/purecn/log/panel_unknown")
    actual = panel_of_normals_workflow.get_log_file("purecn", "prepare")
    assert actual == expected


def test_purecn_step_part_get_input_files_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "container": "work/containers/out/purecn.simg",
        "intervals": "work/purecn/out/panel_unknown.list",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
    }
    actual = panel_of_normals_workflow.get_input_files("purecn", "coverage")(wildcards)
    assert actual == expected


def test_purecn_step_part_get_output_files_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_coverage()"""
    expected = {
        "coverage": "work/{mapper}.purecn/out/{mapper}.purecn.{library_name,.+-DNA[0-9]+-WES[0-9]+}_coverage_loess.txt.gz"
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_log_file_coverage(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_coverage()"""
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.purecn/log/{mapper}.purecn.{library_name,.+-DNA[0-9]+-WES[0-9]+}"
    )
    actual = panel_of_normals_workflow.get_log_file("purecn", "coverage")
    assert actual == expected


def test_purecn_step_part_get_input_files_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "normals": [
            "work/bwa.purecn/out/bwa.purecn.P001-N1-DNA1-WGS1_coverage_loess.txt.gz",
            "work/bwa.purecn/out/bwa.purecn.P002-N1-DNA1-WGS1_coverage_loess.txt.gz",
        ],
        "container": "work/containers/out/purecn.simg",
    }
    actual = panel_of_normals_workflow.get_input_files("purecn", "create_panel")(wildcards)
    assert actual == expected


def test_purecn_step_part_get_output_files_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_output_files_create_panel()"""
    expected = {
        "db": "work/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.rds",
        "db_md5": "work/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.rds.md5",
        "mapbias": "work/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.rds",
        "mapbias_md5": "work/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.rds.md5",
        "lowcov": "work/{mapper}.purecn/out/{mapper}.purecn.low_coverage_targets.bed",
        "hq": "work/{mapper}.purecn/out/{mapper}.purecn.hq_sites.bed",
        "plot": "work/{mapper}.purecn/out/{mapper}.purecn.interval_weights.png",
    }
    actual = panel_of_normals_workflow.get_output_files("purecn", "create_panel")
    assert actual == expected


def test_purecn_step_part_get_log_file_create_panel(panel_of_normals_workflow):
    """Tests PureCnStepPart._get_log_file_create_panel()"""
    expected = get_expected_log_files_dict(
        base_out="work/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals"
    )
    actual = panel_of_normals_workflow.get_log_file("purecn", "create_panel")
    assert actual == expected


def test_purecn_step_part_get_resource_usage(panel_of_normals_workflow):
    """Tests PureCnStepPart.get_resource_usage() for all actions"""
    expected = {
        "coverage": {"threads": 1, "memory": "24G", "time": "04:00:00"},
        "prepare": {"threads": 1, "memory": "24G", "time": "04:00:00"},
        "create_panel": {"threads": 1, "memory": "32G", "time": "12:00:00"},
    }
    for action, resources in expected.items():
        for resource, value in resources.items():
            actual = panel_of_normals_workflow.get_resource("purecn", action, resource)()
            assert actual == value


# Tests for CnvkitStepPart ------------------------------------------------------------------------


def test_cnvkit_step_part_get_args_access(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_args_access()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    actual = panel_of_normals_workflow.get_args("cnvkit", "access")(
        wildcards,
        panel_of_normals_workflow.get_input_files("cnvkit", "access")(wildcards),
    )
    if actual.get("ignore_chroms", None) is not None:
        actual["ignore_chroms"].sort()
    expected = {"reference": "/path/to/ref.fa", "min-gap-size": None, "exclude": [], "ignore_chroms": ["GL*", "MT"]}
    assert actual == expected


def test_cnvkit_step_part_get_args_autobin(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_args_access()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "bams": [
            "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
            "NGS_MAPPING/output/bwa.P002-N1-DNA1-WGS1/out/bwa.P002-N1-DNA1-WGS1.bam",
        ],
        "access": "work/bwa.cnvkit/out/bwa.cnvkit.access.bed",
        "method": "wgs",
        "bp-per-bin": 50000,
    }
    actual = panel_of_normals_workflow.get_args("cnvkit", "autobin")(
        wildcards,
        panel_of_normals_workflow.get_input_files("cnvkit", "autobin")(wildcards),
    )
    assert actual == expected


def test_cnvkit_step_part_get_args_target(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_args_target()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "interval": "work/bwa.cnvkit/out/bwa.cnvkit.access.bed",
        "avg-size": 2000,
        "split": True,
        "annotate": "/path/to/annotations.gtf",
        "short-names": True,

    }
    actual = panel_of_normals_workflow.get_args("cnvkit", "target")(
        wildcards,
        panel_of_normals_workflow.get_input_files("cnvkit", "target")(wildcards),
    )
    assert actual == expected


def test_cnvkit_step_part_get_args_coverage(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_args_coverage()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "region": "target",
        }
    )
    expected = {
        "intervals": "work/bwa.cnvkit/out/bwa.cnvkit.target.bed",
        "bam": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam.bai",
        "reference": "/path/to/ref.fa",
        "min-mapq": 0,
        "count": False,
    }
    actual = panel_of_normals_workflow.get_args("cnvkit", "coverage")(
        wildcards,
        panel_of_normals_workflow.get_input_files("cnvkit", "coverage")(wildcards),
    )
    assert actual == expected


def test_cnvkit_step_part_get_args_reference(panel_of_normals_workflow):
    """Tests CnvkitStepPart._get_args_create_panel()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    expected = {
        "normals": [
            "work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.target.cnn",
            "work/bwa.cnvkit.P002-N1-DNA1-WGS1/out/bwa.cnvkit.P002-N1-DNA1-WGS1.target.cnn",
        ],
        "reference": "/path/to/ref.fa",
        "cluster": False,
        "no-gc": False,
        "no-edge": True,
        "no-rmask": False,
        "male-reference": False,
        "diploid-parx-genome": "GRCh38",
    }
    actual = panel_of_normals_workflow.get_args("cnvkit", "create_panel")(
        wildcards,
        panel_of_normals_workflow.get_input_files("cnvkit", "create_panel")(wildcards),
    )
    assert actual == expected


def test_cnvkit_step_parts_get_output_files(panel_of_normals_workflow):
    """Tests CnvkitStepPart.get_output_files() for all actions"""
    actions = {
        "access": {"access": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.access.bed"},
        "autobin": {"result": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.autobin.txt"},
        "target": {"target": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.target.bed"},
        "antitarget": {"antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed"},
        "coverage": {"coverage": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{region,(target|antitarget)}.cnn"},
        "create_panel": {"reference": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.panel_of_normals.cnn"},
        "sex": {"sex": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.sex.tsv"},
    }
    for action, result in actions.items():
        expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
        actual = panel_of_normals_workflow.get_output_files("cnvkit", action)
        assert actual == expected


def test_cnvkit_step_parts_get_log_file(panel_of_normals_workflow):
    """Tests CnvkitStepPart.get_log_file() for all actions"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    actions = ("autobin", "target", "create_panel", "sex")
    base_log = "work/{mapper}.cnvkit/log/{mapper}.cnvkit"
    for action in actions:
        result = {k: base_log + f".{action}.{v}" for k, v in exts} 
        expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
        actual = panel_of_normals_workflow.get_log_file("cnvkit", action)
        assert actual == expected


def test_cnvkit_step_parts_get_log_file_access(panel_of_normals_workflow):
    """Tests CnvkitStepPart.get_log_file() for access"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    base_log = "work/{mapper}.cnvkit/log/{mapper}.cnvkit.access"
    result = {k: base_log + f".{v}" for k, v in exts} 
    expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "access")
    assert actual == expected


def test_cnvkit_step_parts_get_log_file_coverage(panel_of_normals_workflow):
    """Tests CnvkitStepPart.get_log_file() for coverage"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    base_log = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}.{region,(target|antitarget)}"
    result = {k: base_log + f".{v}" for k, v in exts} 
    expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
    actual = panel_of_normals_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


# PanelOfNormalsWorkflow  --------------------------------------------------------------------------


def test_panel_of_normals_workflow(panel_of_normals_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvkit", "link_out", "mutect2", "purecn"]
    actual = list(sorted(panel_of_normals_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    expected = []
    tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.panel_of_normals.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext)
        for ext in ("vcf.gz", "vcf.gz.md5", "vcf.gz.tbi", "vcf.gz.tbi.md5")
        for mapper in ("bwa",)
    ]
    tpl = "output/{mapper}.mutect2/out/{mapper}.mutect2.genomicsDB.{ext}"
    expected += [
        tpl.format(mapper=mapper, ext=ext)
        for ext in ("tar.gz", "tar.gz.md5")
        for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.mutect2/log/{mapper}.mutect2.panel_of_normals"
    for mapper in ("bwa",):
        expected += get_expected_log_files_dict(base_out=tpl.format(mapper=mapper)).values()

    # Now for basic cnvkit files (panel of normal only)
    tpl = "output/{mapper}.cnvkit/out/{mapper}.cnvkit.{substep}.{ext}{chksum}"
    expected += [
        tpl.format(mapper=mapper, substep=substep, ext=ext, chksum=chksum)
        for chksum in ("", ".md5")
        for (substep, ext) in (("panel_of_normals", "cnn"), ("sex", "tsv"), ("target", "bed"))
        for mapper in ("bwa",)
    ]
    # add log files
    tpl = "output/{mapper}.cnvkit/log/{mapper}.cnvkit.{substep}"
    for substep in ("create_panel", "sex", "target"):
        for mapper in ("bwa",):
            base_out = tpl.format(mapper=mapper, substep=substep)
            expected += get_expected_log_files_dict(base_out=base_out).values()
            expected += [base_out + ".sh", base_out + ".sh.md5"]

    # PureCN
    tpl = "output/{mapper}.purecn/out/{mapper}.purecn.panel_of_normals.rds{chksum}"
    expected += [
        tpl.format(mapper=mapper, chksum=chksum) for mapper in ("bwa",) for chksum in ("", ".md5")
    ]
    tpl = "output/{mapper}.purecn/out/{mapper}.purecn.mapping_bias.rds{chksum}"
    expected += [
        tpl.format(mapper=mapper, chksum=chksum) for mapper in ("bwa",) for chksum in ("", ".md5")
    ]
    expected += get_expected_log_files_dict(
        base_out="output/{mapper}.purecn/log/{mapper}.purecn.panel_of_normals".format(mapper="bwa")
    ).values()
    tpl = "output/purecn/out/panel_unknown.{ext}{chksum}"
    expected += [
        tpl.format(ext=ext, chksum=chksum)
        for ext in ("list", "bed.gz", "bed.gz.tbi")
        for chksum in ("", ".md5")
    ]
    expected += get_expected_log_files_dict(
        base_out="output/purecn/log/panel_unknown".format()
    ).values()

    expected = list(sorted(expected))
    actual = list(sorted(panel_of_normals_workflow.get_result_files()))
    assert actual == expected

# -*- coding: utf-8 -*-
"""Tests for the panel_of_normals workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_cnv_calling import SomaticCnvCallingWorkflow

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

          somatic_variant_calling:
            tools: ['mutect2']
            path_ngs_mapping: ../ngs_mapping
            mutect2:
              common_variants: /path/to/common/variants
        
          somatic_purity_ploidy_estimate:
            tools: ['ascat']
            path_ngs_mapping: ../ngs_mapping
            ascat:
              b_af_loci: /path/to/locii.bed

          somatic_cnv_calling:
              ignore_chroms: [GL*]
              tools:
                wgs: ['cnvkit']
              path_ngs_mapping: ../ngs_mapping
              cnvkit:
                ignore_chroms: [MT]
                diploid_parx_genome: GRCh38
                panel_of_normals:
                  source: paired
                somatic_variant_calling:
                  enabled: True
                  source: cohort
                  tool: mutect2
                  path_somatic_variant_calling: ../somatic_variant_calling
                somatic_purity_ploidy_estimate:
                  enabled: True
                  source: cohort
                  tool: ascat
                segment:
                  threshold: 0.0001
                scatter:
                  enabled: true

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
def somatic_cnv_calling_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    autobin_result_calling_fake_fs,
    purity_result_fake_fs,
    aligner_indices_fake_fs,
    mocker,
):
    """Return PanelOfNormalsWorkflow object pre-configured with germline sheet"""
    # Patch out file-system to enable reading autobin output
    autobin_result_calling_fake_fs.fs.create_file(
        file_path="work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.autobin.txt",
        contents="Target: -1 2000\n",
        create_missing_dirs=True,
    )
    # Patch out file-system to enable reading autobin output
    purity_result_fake_fs.fs.create_file(
        file_path="SOMATIC_PURITY_PLOIDY_ESTIMATE/output/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.txt",
        contents="Purity/ploidy:\t0.35\t2.2\n",
        create_missing_dirs=True,
    )
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.somatic_cnv_calling", autobin_result_calling_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there
    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_calling_cnvkit": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
        "panel_of_normals_cnvkit": lambda x: "SOMATIC_VARIANT_CALLING/" + x,
        "somatic_purity_ploidy_estimate_cnvkit": lambda x: "SOMATIC_PURITY_PLOIDY_ESTIMATE/" + x,
    }
    # Construct the workflow object
    return SomaticCnvCallingWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


# Tests for CnvkitStepPart ------------------------------------------------------------------------


def test_cnvkit_step_part_get_input_files_ignore(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_ignore()"""
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "ignore")(None)
    expected = {"reference": "/path/to/ref.fa"}
    assert actual == expected


def test_cnvkit_step_part_get_args_ignore(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_ignore()"""
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "ignore")(None, None)
    actual["ignore_chroms"].sort()
    expected = {"ignore_chroms": ["GL*", "MT"]}
    assert actual == expected


def test_cnvkit_step_part_get_input_files_access(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_access()"""
    wildcards = Wildcards(fromdict={"mapper": "bwa"})
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "access")(wildcards)
    expected = {
        "reference": "/path/to/ref.fa",
        "exclude": [],
        "ignore_chroms": "work/bwa.cnvkit/out/bwa.cnvkit.ignored.bed",
    }
    assert actual == expected


def test_cnvkit_step_part_get_args_access(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_access()"""
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "access")(None, None)
    expected = {"min-gap-size": None}
    assert actual == expected


def test_cnvkit_step_part_get_input_files_autobin(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_access()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "bams": ["NGS_MAPPING/output/bwa.P001-N1-DNA1-WGS1/out/bwa.P001-N1-DNA1-WGS1.bam"],
        "access": "work/bwa.cnvkit/out/bwa.cnvkit.access.bed",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "autobin")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_autobin(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_access()"""
    expected = {"method": "wgs", "bp-per-bin": 50000}
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "autobin")(None, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_target(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_target()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {
        "interval": "work/bwa.cnvkit/out/bwa.cnvkit.access.bed",
        "avg_size": "work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.autobin.txt",
        "annotate": "/path/to/annotations.gtf",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "target")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_target(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_target()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
        }
    )
    expected = {"avg-size": 2000, "split": True, "short-names": True}
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "target")(
        wildcards,
        somatic_cnv_calling_workflow.get_input_files("cnvkit", "target")(wildcards)
    )
    assert actual == expected


def test_cnvkit_step_part_get_input_files_coverage(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_coverage()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
            "region": "target",
        }
    )
    expected = {
        "intervals": "work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.target.bed",
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "bai": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam.bai",
        "reference": "/path/to/ref.fa",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "coverage")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_coverage(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_coverage()"""
    expected = {"min-mapq": 0, "count": False}
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "coverage")(None, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_reference(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-N1-DNA1-WGS1",
            "region": "target",
        }
    )
    expected = {
        "normals": [
            "work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.targetcoverage.cnn",
        ],
        "reference": "/path/to/ref.fa",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "reference")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_reference(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_create_panel()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "cluster": False,
        "no-gc": False,
        "no-edge": True,
        "no-rmask": False,
        "male-reference": False,
        "diploid-parx-genome": "GRCh38",
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "reference")(wildcards, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_fix(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_fix()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "target": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.targetcoverage.cnn",
        "reference": "work/bwa.cnvkit.P001-N1-DNA1-WGS1/out/bwa.cnvkit.P001-N1-DNA1-WGS1.reference.cnn",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "fix")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_fix(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_fix()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "cluster": False,
        "no-gc": False,
        "no-edge": True,
        "no-rmask": False,
        "diploid-parx-genome": "GRCh38",
        "sample-id": "P001-T1-DNA1-WGS1",
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "fix")(wildcards, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_segment(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_segment()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "ratios": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr",
        "variants": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.vcf.gz",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "segment")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_segment(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_segment()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "method": "cbs",
        "threshold": 0.0001,
        "smooth-cbs": False,
        "drop-low-coverage": False,
        "drop-outliers": 10,
        "sample-id": "P001-T1-DNA1-WGS1",
        "normal-id": "P001-N1-DNA1-WGS1",
        "min-variant-depth": 20,
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "segment")(wildcards, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_call(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_call()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "segments": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/report/bwa.cnvkit.P001-T1-DNA1-WGS1.segmetrics.cns",
        "variants": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.vcf.gz",
        "purity_file": "SOMATIC_PURITY_PLOIDY_ESTIMATE/output/bwa.ascat.P001-T1-DNA1-WGS1/out/bwa.ascat.P001-T1-DNA1-WGS1.txt",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "call")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_call(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_call()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "method": "threshold",
        "thresholds": [-1.1, -0.25, 0.2, 0.7],
        "drop-low-coverage": False,
        "male-reference": False,
        "diploid-parx-genome": "GRCh38",
        "sample-id": "P001-T1-DNA1-WGS1",
        "normal-id": "P001-N1-DNA1-WGS1",
        "min-variant-depth": 20,
        "purity": 0.35,
        "ploidy": 2.2,
        "filter": ["ci"],
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "call")(
        wildcards,
        somatic_cnv_calling_workflow.get_input_files("cnvkit", "call")(wildcards),
    )
    assert actual == expected


def test_cnvkit_step_part_get_input_files_bintest(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_bintest()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "segments": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns",
        "ratios": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "bintest")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_bintest(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_bintest()"""
    expected = {
        "alpha": 0.005,
        "target": False,
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "bintest")(None, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_metrics(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_metrics()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "segments": ["work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns"],
        "ratios": ["work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr"],
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "metrics")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_metrics(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_metrics()"""
    expected = {"drop-low-coverage": False}
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "metrics")(None, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_segmetrics(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_segmetrics()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "segments": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns",
        "ratios": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "segmetrics")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_segmetrics(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_segmetrics()"""
    expected = {
        "drop-low-coverage": False,
        "alpha": 0.05,
        "bootstrap": 100,
        "smooth-bootstrap": True,
        "stats": ["mean", "median", "mode", "t-test", "stdev", "sem", "mad", "mse", "iqr", "bivar", "ci", "pi"]
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "segmetrics")(None, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_genemetric(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_genemetrics()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "segments": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns",
        "ratios": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "genemetrics")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_genemetric(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_genemetrics()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "threshold": 0.2,
        "min-probes": 3,
        "drop-low-coverage": False,
        "male-reference": False,
        "diploid-parx-genome": "GRCh38",
        "alpha": 0.05,
        "bootstrap": 100,
        "stats": ["mean", "median", "mode", "ttest", "stdev", "sem", "mad", "mse", "iqr", "bivar", "ci", "pi"]
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "genemetrics")(wildcards, None)
    assert actual == expected


def test_cnvkit_step_part_get_input_files_scatter(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_input_files_scatter()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
            "contig_name": "1",
        }
    )
    expected = {
        "segments": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cns",
        "ratios": "work/bwa.cnvkit.P001-T1-DNA1-WGS1/out/bwa.cnvkit.P001-T1-DNA1-WGS1.cnr",
        "variants": "SOMATIC_VARIANT_CALLING/output/bwa.mutect2.P001-T1-DNA1-WGS1/out/bwa.mutect2.P001-T1-DNA1-WGS1.vcf.gz",
    }
    actual = somatic_cnv_calling_workflow.get_input_files("cnvkit", "scatter")(wildcards)
    assert actual == expected


def test_cnvkit_step_part_get_args_scatter(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart._get_args_scatter()"""
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "library_name": "P001-T1-DNA1-WGS1",
            "contig_name": "1",
        }
    )
    expected = {
        "chromosome": "1",
        "width": 1000000,
        "antitarget-marker": "o",
        "by-bin": False,
        "trend": False,
        "segment-color": "darkorange",
        "title": "P001-T1-DNA1-WGS1 - 1",
        "sample-id": "P001-T1-DNA1-WGS1",
        "normal-id": "P001-N1-DNA1-WGS1",
        "min-variant-depth": 20,
        "fig-size": (12.256, 16.192),
    }
    actual = somatic_cnv_calling_workflow.get_args("cnvkit", "scatter")(wildcards, None)
    assert actual == expected


def test_cnvkit_step_parts_get_output_files(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart.get_output_files() for all actions"""
    actions = {
        "ignore": {"ignore_chroms": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.ignored.bed"},
        "access": {"access": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.access.bed"},
        "autobin": {"result": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.autobin.txt"},
        "target": {"target": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.target.bed"},
        "antitarget": {"antitarget": "work/{mapper}.cnvkit/out/{mapper}.cnvkit.antitarget.bed"},
        "coverage": {"coverage": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{region,(target|antitarget)}coverage.cnn"},
        "reference": {"reference": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.reference.cnn"},
        "fix": {"ratios": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cnr"},
        "segment": {
            "segments": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.cns",
            "dataframe": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.rds",
        },
        "call": {"calls": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.call.cns"},
        "bintest": {"tests": "work/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.bintest.cns"},
        "metrics": {"report": "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}.metrics.tsv"},
        "genemetrics": {"report": "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}.genemetrics.tsv"},
        "segmetrics": {"report": "work/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}.segmetrics.cns"},
        "scatter": {"plot": "work/{mapper}.cnvkit.{library_name}/plot/{mapper}.cnvkit.{library_name}.scatter.{contig_name}.jpeg"},
    }
    for action, result in actions.items():
        expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
        actual = somatic_cnv_calling_workflow.get_output_files("cnvkit", action)
        assert actual == expected


def test_cnvkit_step_parts_get_log_file(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart.get_log_file() for all actions"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    actions = ("autobin", "target", "reference", "fix", "segment", "call", "bintest", "metrics", "segmetrics", "genemetrics")
    base_log = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}"
    for action in actions:
        result = {k: base_log + f".{action}.{v}" for k, v in exts} 
        expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
        actual = somatic_cnv_calling_workflow.get_log_file("cnvkit", action)
        assert actual == expected


def test_cnvkit_step_parts_get_log_file_no_lib(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart.get_log_file() for all actions not dependent on library"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    actions = ("ignore", "access")
    base_log = "work/{mapper}.cnvkit/log/{mapper}.cnvkit."
    for action in actions:
        result = {k: base_log + f"{action}.{v}" for k, v in exts} 
        expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
        actual = somatic_cnv_calling_workflow.get_log_file("cnvkit", action)
        assert actual == expected


def test_cnvkit_step_parts_get_log_file_coverage(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart.get_log_file() for coverage"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    base_log = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}.{region,(target|antitarget)}coverage"
    result = {k: base_log + f".{v}" for k, v in exts} 
    expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
    actual = somatic_cnv_calling_workflow.get_log_file("cnvkit", "coverage")
    assert actual == expected


def test_cnvkit_step_parts_get_log_file_scatter(somatic_cnv_calling_workflow):
    """Tests CnvkitStepPart.get_log_file() for coverage"""
    exts = (("conda_info", "conda_info.txt"), ("conda_list", "conda_list.txt"), ("log", "log"), ("sh", "sh"))
    base_log = "work/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}.scatter.{contig_name}"
    result = {k: base_log + f".{v}" for k, v in exts} 
    expected = result | {k + "_md5": v + ".md5" for k, v in result.items()}
    actual = somatic_cnv_calling_workflow.get_log_file("cnvkit", "scatter")
    assert actual == expected


# SomaticCnvCallingWorkflow  --------------------------------------------------------------------------


def test_somatic_cnv_calling_workflow(somatic_cnv_calling_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["cnvkit", "link_out"]
    actual = list(sorted(somatic_cnv_calling_workflow.sub_steps.keys()))
    assert actual == expected

    tumor_libraries = ("P001-T1-DNA1-WGS1", "P002-T1-DNA1-WGS1", "P002-T2-DNA1-WGS1")

    expected = []

    # cnvkit output files
    tpl = "output/{mapper}.cnvkit.{library_name}/out/{mapper}.cnvkit.{library_name}.{ext}"
    expected += [
        tpl.format(mapper=mapper, library_name=library_name, ext=ext)
        for ext in ("cnr", "cns", "call.cns", "bintest.cns")
        for library_name in tumor_libraries
        for mapper in ("bwa",)
    ]

    # cnvkit log files
    tpl = "output/{mapper}.cnvkit.{library_name}/log/{mapper}.cnvkit.{library_name}.{step}.{ext}"
    expected += [
        tpl.format(mapper=mapper, library_name=library_name, step=step, ext=ext)
        for ext in ("conda_info.txt", "conda_list.txt", "log", "sh")
        for step in ("fix", "segment", "call", "bintest")
        for library_name in tumor_libraries
        for mapper in ("bwa",)
    ]

    # cnvkit report files
    tpl = "output/{mapper}.cnvkit.{library_name}/report/{mapper}.cnvkit.{library_name}.{step}.{ext}"
    expected += [
        tpl.format(mapper=mapper, library_name=library_name, step=step, ext=ext)
        for step, ext in (("metrics", "tsv"), ("genemetrics", "tsv"), ("segmetrics", "cns"))
        for library_name in tumor_libraries
        for mapper in ("bwa",)
    ]

    # cnvkit plot files
    tpl = "output/{mapper}.cnvkit.{library_name}/plot/{mapper}.cnvkit.{library_name}.{step}.{contig_name}.{ext}"
    expected += [
        tpl.format(mapper=mapper, library_name=library_name, step=step, contig_name=contig_name, ext=ext)
        for ext in ("jpeg",)
        for contig_name in ["all"] + list(map(str, range(1, 23))) + ["X", "Y"]
        for step in ("scatter",)
        for library_name in tumor_libraries
        for mapper in ("bwa",)
    ]

    # Add md5
    expected += [x + ".md5" for x in expected]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_cnv_calling_workflow.get_result_files()))
    assert actual == expected

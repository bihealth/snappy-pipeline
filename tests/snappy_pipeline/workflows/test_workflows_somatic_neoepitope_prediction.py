# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_neoepitope_prediction import (
    SomaticNeoepitopePredictionWorkflow,
)

from .common import get_expected_log_files_dict, get_expected_output_vcf_files_dict
from .conftest import patch_module_fs


# Test tumor mutational burden calculation with vcf file from somatic variant calling step
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
          cosmic:
            path: /path/to/cosmic.vcf.gz
          dbsnp:
            path: /path/to/dbsnp.vcf.gz
          features:
            path: /path/to/genecode.gtf

        step_config:
          ngs_mapping:
            tools:
              dna: ['bwa']
              rna: ['star']
            star:
              path_index: /path/to/star/index
            bwa:
              path_index: /path/to/bwa/index.fa

          gene_expression_quantification:
            tools: [salmon]
            salmon:
              path_index: /path/to/salmon/index

          somatic_variant_calling:
            tools: [mutect2]
            mutect2:
              contamination:
                common_variants: /path/to/common/variants

          hla_typing:
            path_link_in: ''  # OPTIONAL Override data set configuration search paths for FASTQ files
            tools: [optitype]   # REQUIRED - available: 'optitype' and 'arcashla'
            optitype:
                max_reads: 5000
                num_mapping_threads: 4

          somatic_variant_annotation:
            path_somatic_variant: ../somatic_variant_calling
            tools: ["vep"]
            vep:
                cache_dir: /path/to/dir/cache

          somatic_neoepitope_prediction:
            tools: [pVACseq]
            pileup:
              enabled: true
            quantification:
              enabled: true
            pVACseq:
                algorithms: ['MHCflurry','MHCnuggetsI']
                
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
def somatic_neoepitope_prediction_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    hla_typing_result_fake_fs,
    mocker,
):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping.model", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.hla_typing.model", aligner_indices_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there

    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_annotation": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
        "hla_typing": lambda x: "HLA_TYPING/" + x,
        "gene_expression_quantification": lambda x: "GENE_EXPRESSION_QUANTIFICATION/" + x
    }
    # Construct the workflow object
    return SomaticNeoepitopePredictionWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )


def test_somatic_neoepitope_pileup_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    # Define expected
    annotated_tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    mapping_tpl = "{mapper}.P001-T1-RNA1-mRNA_seq1"
    expected = {
        "bam": f"NGS_MAPPING/output/{mapping_tpl}/out/{mapping_tpl}.bam",
        "loci": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.vcf.gz",
        "reference": "/path/to/ref.fa",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "pileup")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_combine_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    # Define expected
    annotated_tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expr_tpl = "{mapper}.P001-T1-RNA1-mRNA_seq1"
    expected = {
        "pileup": f"work/{annotated_tpl}/out/{annotated_tpl}.pileup.vcf.gz",
        "annotated": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.vcf.gz",
        "gene_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.gene.sf",
        "transcript_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.transcript.sf",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "combine")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    annotated_tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    hla_tpl = "HLA_TYPING/output/{{typer}}.{library}/out/{{typer}}.{library}.txt"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "vcf": f"work/{annotated_tpl}/out/{annotated_tpl}.combined.vcf.gz",
        "hla_tumor_dna": hla_tpl.format(library="P001-T1-DNA1-WGS1"),
        "hla_normal_dna": hla_tpl.format(library="P001-N1-DNA1-WGS1"),
        "hla_tumor_rna": hla_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "predict")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_install_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    expected = {"container": "work/containers/out/pvactools.sif"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "install")
    assert actual == expected


def test_somatic_neoepitope_pileup_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.pileup.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "pileup")
    assert actual == expected


def test_somatic_neoepitope_combine_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.combined.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "combine")
    assert actual == expected


def test_somatic_neoepitope_predict_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = {"done": f"work/{tpl}/neoepitopes/.done"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "predict")
    assert actual == expected


def test_somatic_neoepitope_preparation_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {
        "install": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "pileup": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "combine": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "predict": {"threads": 4, "time": "23:59:59", "memory": "64G"},
    }
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected} of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("pvacseq", action, resource)()
            assert actual == expected, msg_error


def test_somatic_neoepitope_install_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    base_out = "work/containers/log/pvactools"
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "install")
    assert actual == expected


def test_somatic_neoepitope_pileup_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/pileup")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "pileup")
    assert actual == expected


def test_somatic_neoepitope_combine_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/combine")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "combine")
    assert actual == expected


def test_somatic_neoepitope_predict_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.filtered.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/predict")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "predict")
    assert actual == expected


def test_somatic_neoepitope_pileup_step_part_get_args(somatic_neoepitope_prediction_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "extra_args": "--adjust-MQ 0 --delta-BQ 30 --max-BQ 60 --max-depth 250 --min-BQ 1 --min-MQ 0",
        "tumor_sample": "P001-T1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "pileup")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_combine_step_part_get_args(somatic_neoepitope_prediction_workflow):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "tumor_sample": "P001-T1",
        "normal_sample": "P001-N1",
        "tumor_library": "P001-T1-DNA1-WGS1",
        "normal_library": "P001-N1-DNA1-WGS1",
        "extra_args": (
            r"--ensembl-id --annotation 'CSQ' --annotation-description-regex "
            r"'^Consequence annotations from Ensembl VEP. Format: (?P<titles>.+)$' "
            r"--annotation-gene-id 'Gene' --annotation-separator '\|' "
            r"--annotation-transcript-id 'Feature' --format 'salmon'"
        ),
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "combine")(wildcards)
    assert actual == expected


# def test_somatic_neoepitope_predict_step_part_get_args(somatic_neoepitope_prediction_workflow):
#     wildcards = Wildcards(
#         fromdict={
#             "mapper": "bwa",
#             "caller": "mutect2",
#             "annotator": "vep",
#             "tumor_dna": "P001-T1-DNA1-WGS1",
#         }
#     )
#     expected = {
#         "num_threads": 1,
#         "tumor_sample": "P001-T1",
#         "normal_sample": "P001-N1",
#         "hla_type_I": ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*15:32", "HLA-C*04:03"],
#         "extra_args": (
#             "--aggregate-inclusion-binding-threshold 5000 --aggregate-inclusion-count-limit 25 "
#             "--algorithms 'all_class_i' --class_i-epitope-length '8,9,10,11' "
#             "--downstream-sequence-length 1000 --expn-val 1.0 --minimum-fold-change 0.0 "
#             "--netmhciipan-version '4.1' --netmhc_stab --net-chop-method 'cterm' "
#             "--normal-cov 25 --normal-vaf 0.02 --precentile-threshold 1.0 "
#             "--percentage-threshold-strategy 'conservative' --tdna-cov 25 --tdna-vaf 0.1 "
#             "--top-score-metric 'median' --top-score-metric2 'percentile' "
#             "--transcript-prioritization-strategy 'mane_select' --trna-cov 2 --trna-vaf 0.25"
#         )
#     }
#     actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "predict")(wildcards)
#     assert actual == expected


def test_somatic_neoepitope_prediction_workflow(somatic_neoepitope_prediction_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "pvacseq"]
    actual = list(sorted(somatic_neoepitope_prediction_workflow.sub_steps.keys()))
    assert actual == expected

    # Check result file construction
    tpl = "output/{mapper}.{caller}.{annotator}.filtered.P00{i}-T{t}-DNA1-WGS1/log/predict.log"
    expected = [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", i=i, t=t)
        for (i, t) in ((1, 1), (2,1), (2, 2))
    ]

    tpl = "output/{mapper}.{caller}.{annotator}.filtered.P00{i}-T{t}-DNA1-WGS1/neoepitopes/.done"
    expected += [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", i=i, t=t)
        for (i, t) in ((1, 1), (2, 1), (2, 2))
    ]
    expected = list(sorted(expected))
    actual = list(sorted(somatic_neoepitope_prediction_workflow.get_result_files()))
    assert expected == actual

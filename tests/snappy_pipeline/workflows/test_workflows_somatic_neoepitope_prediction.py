# -*- coding: utf-8 -*-
"""Tests for the somatic_variant_calling workflow module code"""

import textwrap

import pytest
import ruamel.yaml as ruamel_yaml
from snakemake.io import Wildcards

from snappy_pipeline.workflows.somatic_neoepitope_prediction import (
    SomaticNeoepitopePredictionWorkflow,
)

from .common import get_expected_log_files_dict
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
            path: /path/to/gencode.gtf

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
            tools:
              dna: [optitype]   # REQUIRED - available: 'optitype' and 'arcashla'
              rna: [optitype, arcashla]
            optitype:
                max_reads: 5000
                num_mapping_threads: 4

          somatic_variant_annotation:
            path_somatic_variant: ../somatic_variant_calling
            tools: ["vep"]
            is_filtered: false
            vep:
                cache_dir: /path/to/dir/cache

          germline_variant_calling:
            path_ngs_mapping: NGS_MAPPING
            tools: [gatk4_hc]
            gatk4_hc:
              num_threads: 8

          germline_variant_filtration:
            path_variant: ../germline_variant_calling
            has_annotation: false
            filter_list:
              - bcftools:
                  exclude: "AD[0:0]+AD[0:1]<50 | AD[0:1]<5 | AD[0:1]/(AD[0:0]+AD[0:1])<0.05"

          combine_variants:
            somatic_variant_type: annotation
            path_somatic_variant: SOMATIC_VARIANT_ANNOTATION
            tool_somatic_variant_annotation: vep
            germline_variant_type: filtration
            path_germline_variant: GERMLINE_VARIANT_FILTRATION
            is_germline_variant_filtered: true
            rename_combined: tumor

          somatic_neoepitope_prediction:
            tools: [pvacseq, pvacfuse, pvacsplice]
            is_filtered: false
            pileup:
              enabled: true
            quantification:
              enabled: true
              duplicate_transcripts_table: /path/to/duplicates.tsv
            phasing:
              enabled: true
              path_combine_variants: COMBINE_VARIANTS
            proteome:
              enabled: true
              path_germline_variants: GERMLINE_VARIANT_FILTRATION
              germline_variant_step: germline_variant_filtration
              external_proteome: /path/to/gencode.fa
            tools_hla_typing:
              dna:
                class_i: [optitype]
              rna:
                class_i: [optitype]
                class_ii: [arcashla]
            pvacseq:
                class_ii_epitope_length: [10, 11]
                extra_args: ["--percentile-threshold-strategy exploratory"]
                net_chop:
                  enabled: true
                  path_netchop: /path/to/netchop.bin
                netmhc_stab:
                  enabled: true
            pvacfuse:
                algorithms: [all_class_i]
            pvacsplice:
                algorithms: [NetMHCIIpan, MHCnuggetsII]
                use_all_transcripts: true
                genes_of_interest_file: /path/to/genes.txt
                net_chop:
                  enabled: true
                  path_netchop: /path/to/netchop.bin
                netmhc_stab:
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
def somatic_neoepitope_prediction_workflow(
    dummy_workflow,
    minimal_config,
    config_lookup_paths,
    work_dir,
    config_paths,
    cancer_sheet_fake_fs,
    aligner_indices_fake_fs,
    hla_typing_result_fake_fs,
    strandedness_result_fake_fs,
    mocker,
):
    # Patch out file-system related things in abstract (the crawling link in step is defined there)
    patch_module_fs("snappy_pipeline.workflows.abstract", cancer_sheet_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.ngs_mapping.model", aligner_indices_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.somatic_neoepitope_prediction", hla_typing_result_fake_fs, mocker)
    patch_module_fs("snappy_pipeline.workflows.somatic_neoepitope_prediction", strandedness_result_fake_fs, mocker)
    # Update the "globals" attribute of the mock workflow (snakemake.workflow.Workflow) so we
    # can obtain paths from the function as if we really had a NGSMappingPipelineStep there

    dummy_workflow.globals = {
        "ngs_mapping": lambda x: "NGS_MAPPING/" + x,
        "somatic_variant_annotation": lambda x: "SOMATIC_VARIANT_ANNOTATION/" + x,
        "hla_typing": lambda x: "/HLA_TYPING/" + x,
        "gene_expression_quantification": lambda x: "GENE_EXPRESSION_QUANTIFICATION/" + x,
        "somatic_gene_fusion_calling": lambda x: "SOMATIC_GENE_FUSION_CALLING/" + x,
        "combine_variants": lambda x: "COMBINE_VARIANTS/" + x,
        "germline_variant": lambda x: "GERMLINE_VARIANT_FILTRATION/" + x,
    }
    # Construct the workflow object
    return SomaticNeoepitopePredictionWorkflow(
        dummy_workflow,
        minimal_config,
        config_lookup_paths,
        config_paths,
        work_dir,
    )

# ---- pVACtools (container installation & pvacseq/pvacfuse/pvacsplice output & log)

def test_somatic_neoepitope_prediction_pvactools_install_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    expected = {"container": "work/containers/out/pvactools.sif"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvactools", "install")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_install_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    expected = "work/containers/log/pvactools.log"
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvactools", "install")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_step_part_get_input_files(
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "annotated": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.vcf.gz"
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvactools", "normalize")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.normalized.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvactools", "normalize")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "tumor_library": "P001-T1-DNA1-WGS1",
        "normal_library": "P001-N1-DNA1-WGS1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvactools", "normalize")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/normalize.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvactools", "normalize")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_full_step_part_get_input_files(
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "annotated": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.full.vcf.gz"
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvactools", "normalize_full")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_full_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.normalized.full.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvactools", "normalize_full")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_full_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "tumor_library": "P001-T1-DNA1-WGS1",
        "normal_library": "P001-N1-DNA1-WGS1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvactools", "normalize_full")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvactools_normalize_full_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/normalize_full.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvactools", "normalize_full")
    assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tools = ("pvacseq", "pvacfuse", "pvacsplice")
    exts = (
        ("all", ".all_epitopes.tsv"),
        ("filtered", ".filtered.tsv"),
        ("aggregated", ".all_epitopes.aggregated.tsv"),
        ("json", ".all_epitopes.aggregated.metrics.json"),
    )
    subdirs = (("MHC_Class_I", "MHC_I"), ("MHC_Class_II", "MHC_II"), ("combined", "Combined"))

    for tool in tools:
        tpl = "work/{mapper}.{caller}.{annotator}." + tool + ".{tumor_dna}/out/"
        expected = {"done": tpl + ".done"}
        for k, ext in exts:
            if tool != "pvacseq" and k == "json":
                continue
            for d, fn in subdirs:
                expected[f"{k}.{fn}"] = tpl + f"{d}/{{tumor_dna}}.{fn}{ext}"

        actual = somatic_neoepitope_prediction_workflow.get_output_files(tool, tool)
        assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    for tool in ("pvacseq", "pvacfuse", "pvacsplice"):
        expected = get_expected_log_files_dict(
            base_out="work/{mapper}.{caller}.{annotator}." + tool + ".{tumor_dna}/log/" + tool + ".{tumor_dna}"
        )
        actual = somatic_neoepitope_prediction_workflow.get_log_file(tool, tool)
        assert actual == expected


def test_somatic_neoepitope_proteome_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(fromdict={"tumor_dna": "P001-T1-DNA1-WGS1", "mapper": "bwa", "caller": "gatk4_hc"})
    expected = {
        "vcf": "GERMLINE_VARIANT_FILTRATION/output/bwa.gatk4_hc.filtered.P001-N1-DNA1-WGS1/out/bwa.gatk4_hc.filtered.P001-N1-DNA1-WGS1.vcf.gz",
        "reference": "/path/to/ref.fa",
        "features": "/path/to/gencode.gtf",
        "path_proteome": "/path/to/gencode.fa",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("proteome", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_proteome_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"proteome": f"work/{tpl}/out/{tpl}.proteome.fa.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("proteome", "run")
    assert actual == expected


def test_somatic_neoepitope_proteome_step_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    base_out = f"work/{tpl}/log/proteome.{{tumor_dna}}"
    expected = get_expected_log_files_dict(base_out=base_out)
    actual = somatic_neoepitope_prediction_workflow.get_log_file("proteome", "run")
    assert actual == expected


def test_somatic_neoepitope_proteome_step_part_get_args(
    somatic_neoepitope_prediction_workflow,
):
    expected = {"add_unmutated": True}
    actual = somatic_neoepitope_prediction_workflow.get_args("proteome", "run")
    assert actual == expected


def test_somatic_neoepitope_proteome_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    expected = {"run": {"threads": 1, "time": "03:59:59", "memory": "24G", "partition": "medium"}}
    for action in ("run",):
        for resource in ("threads", "time", "memory", "partition"):
            actual = somatic_neoepitope_prediction_workflow.get_resource("proteome", action, resource)()
            assert actual == expected[action][resource]


# ---- pVACseq

def test_somatic_neoepitope_prediction_pvacseq_pileup_step_part_get_input_files(
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    mapping_tpl = "star.P001-T1-RNA1-mRNA_seq1"
    expected = {
        "bam": f"NGS_MAPPING/output/{mapping_tpl}/out/{mapping_tpl}.bam",
        "loci": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.vcf.gz",
        "reference": "/path/to/ref.fa",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "pileup")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_add_expression_step_part_get_input_files(
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expr_tpl = "salmon.P001-T1-RNA1-mRNA_seq1"
    expected = {
        "pileup": f"work/{annotated_tpl}/out/{annotated_tpl}.pileup.vcf.gz",
        "annotated": f"work/{annotated_tpl}/out/{annotated_tpl}.normalized.vcf.gz",
        "gene_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.gene.sf",
        "transcript_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.transcript.sf",
        "duplicate_transcripts_table": "/path/to/duplicates.tsv",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "add_expression")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_pvacseq_step_part_get_input_files(
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "alleles": f"work/{annotated_tpl}/out/{annotated_tpl}.hla_types.txt",
        "vcf": f"work/{annotated_tpl}/out/{annotated_tpl}.with_expression.vcf.gz",
        "phased": f"work/{annotated_tpl}/out/{annotated_tpl}.phased.vcf.gz",
        "peptides": f"work/{annotated_tpl}/out/{annotated_tpl}.proteome.fa.gz",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "pvacseq")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_pileup_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.pileup.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "pileup")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_add_expression_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.with_expression.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "add_expression")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {
        "pileup": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "add_expression": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "pvacseq": {"threads": 1, "time": "23:59:59", "memory": "64G"},
    }
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected} of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("pvacseq", action, resource)()
            assert actual == expected, msg_error


def test_somatic_neoepitope_prediction_pvacseq_pileup_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/pileup.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "pileup")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_add_expression_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/add_expression.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "add_expression")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_pileup_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "extra_args": [],
        "tumor_sample": "P001-T1-DNA1-WGS1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "pileup")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_add_expression_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "format": "salmon",
        "extra_args": [],
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "add_expression")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_pvacseq_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "algorithms": ["all_class_i", "all_class_ii"],
        "lengths": {
            "class_i": [8, 9, 10, 11],
            "class_ii": [10, 11],
        },
        "exclude_bind": [
            "container",
            "alleles",
            "all.MHC_I",
            "filtered.MHC_I",
            "aggregated.MHC_I",
            "json.MHC_I",
            "all.MHC_II",
            "filtered.MHC_II",
            "aggregated.MHC_II",
            "json.MHC_II",
            "all.Combined",
            "filtered.Combined",
            "aggregated.Combined",
            "json.Combined",
        ],
        "extra_args": [
            "--percentile-threshold-strategy exploratory",
            "--run-ml-predictions",
            "--ml-threshold-accept 0.55",
            "--ml-threshold-reject 0.3",
        ],
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "pvacseq")(wildcards)
    assert actual == expected

# ---- pVACfuse

def test_somatic_neoepitope_prediction_pvacfuse_pvacfuse_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "alleles": f"work/{annotated_tpl}/out/{annotated_tpl}.hla_types.txt",
        "fusions": f"SOMATIC_GENE_FUSION_CALLING/output/arriba.P001-T1-RNA1-mRNA_seq1/out/arriba.P001-T1-RNA1-mRNA_seq1.fusions.tsv",
        "peptides": f"work/{annotated_tpl}/out/{annotated_tpl}.proteome.fa.gz",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacfuse", "pvacfuse")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacfuse_pvacfuse_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "algorithms": ["all_class_i"],
        "lengths": {
            "class_i": [8, 9, 10, 11],
            "class_ii": [12, 13,14,15,16,17, 18],
        },
        "exclude_bind": [
            "container",
            "alleles",
            "all.MHC_I",
            "filtered.MHC_I",
            "aggregated.MHC_I",
            "all.MHC_II",
            "filtered.MHC_II",
            "aggregated.MHC_II",
            "all.Combined",
            "filtered.Combined",
            "aggregated.Combined",
        ],
        "extra_args": [],
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacfuse", "pvacfuse")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacfuse_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"pvacfuse": {"threads": 1, "time": "23:59:59", "memory": "64G"}}
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected} of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("pvacfuse", action, resource)()
            assert actual == expected, msg_error

# ---- pVACsplice

def test_somatic_neoepitope_prediction_pvacsplice_junction_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    rna_tpl = "star.P001-T1-RNA1-mRNA_seq1"
    expected = {
        "annotated": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.full.vcf.gz",
        "bam": f"NGS_MAPPING/output/{rna_tpl}/out/{rna_tpl}.bam",
        "strandedness": f"NGS_MAPPING/output/{rna_tpl}/strandedness/{rna_tpl}.decision.json",
        "reference": "/path/to/ref.fa",
        "features": "work/pvacsplice_workaround/out/features.gtf",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacsplice", "junction")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_pvacsplice_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow
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
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "alleles": f"work/{annotated_tpl}/out/{annotated_tpl}.hla_types.txt",
        "junctions": f"work/{annotated_tpl}/out/{annotated_tpl}.junctions.tsv",
        "annotated": f"work/{annotated_tpl}/out/{annotated_tpl}.normalized.full.vcf.gz",
        "genes": "/path/to/genes.txt",
        "peptides": f"work/{annotated_tpl}/out/{annotated_tpl}.proteome.fa.gz",
        "reference": "/path/to/ref.fa",
        "features": "/path/to/gencode.gtf",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacsplice", "pvacsplice")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_junction_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"junctions": f"work/{tpl}/out/{tpl}.junctions.tsv"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacsplice", "junction")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_junction_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    input = somatic_neoepitope_prediction_workflow.get_input_files("pvacsplice", "junction")(wildcards)
    expected = {"strandedness": "RF"}
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacsplice", "junction")(wildcards, input)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_pvacsplice_step_part_get_args(
    somatic_neoepitope_prediction_workflow
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
        }
    )
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1-DNA1-WGS1",
        "normal_sample": "P001-N1-DNA1-WGS1",
        "lengths": {
            "class_i": [8, 9, 10, 11],
            "class_ii": [12, 13,14,15,16,17, 18],
        },
        "algorithms": ["NetMHCIIpan", "MHCnuggetsII"],
        "exclude_bind": [
            "container",
            "alleles",
            "all.MHC_I",
            "filtered.MHC_I",
            "aggregated.MHC_I",
            "all.MHC_II",
            "filtered.MHC_II",
            "aggregated.MHC_II",
            "all.Combined",
            "filtered.Combined",
            "aggregated.Combined",
        ],
        "extra_args": [],
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacsplice", "pvacsplice")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {
        "junction": {"threads": 1, "time": "03:59:59", "memory": "16G"},
        "pvacsplice": {"threads": 1, "time": "23:59:59", "memory": "64G"},
    }
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected}' of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("pvacsplice", action, resource)()
            assert actual == expected, msg_error

# ---- phasing

def test_somatic_neoepitope_prediction_phasing_step_part_get_input_files(
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
    expected = {
        "bam": "NGS_MAPPING/output/bwa.P001-T1-DNA1-WGS1/out/bwa.P001-T1-DNA1-WGS1.bam",
        "reference": "/path/to/ref.fa",
        "vcf": "COMBINE_VARIANTS/output/bwa.combined.P001-T1-DNA1-WGS1/out/bwa.combined.P001-T1-DNA1-WGS1.vcf.gz",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("phasing", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_phasing_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.phased.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("phasing", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_phasing_step_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/phasing.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("phasing", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_phasing_step_part_get_args(
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
    expected = {}
    actual = somatic_neoepitope_prediction_workflow.get_args("phasing", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_phasing_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"run": {"threads": 1, "time": "23:59:59", "memory": "32G"}}
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected}' of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("phasing", action, resource)()
            assert actual == expected, msg_error

# ---- netchop

def test_somatic_neoepitope_prediction_netchop_step_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
            "tool": "pvacseq",
            "mhc_class_d": "MHC_Class_I",
            "mhc_class_fn": "MHC_I"
        }
    )
    expected = {
        "epitopes": "work/bwa.mutect2.vep.pvacseq.P001-T1-DNA1-WGS1/out/MHC_Class_I/P001-T1-DNA1-WGS1.MHC_I.filtered.tsv",
        "netchop": "/path/to/netchop.bin",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("netchop", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}"
    expected = {"netchop": f"work/{tpl}/out/{{mhc_class_d}}/{{tumor_dna}}.{{mhc_class_fn}}.netchop.tsv"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("netchop", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/netchop.{{mhc_class_d}}_{{mhc_class_fn}}.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("netchop", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_args(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
            "tool": "pvacseq",
            "mhc_class_d": "MHC_Class_I",
            "mhc_class_fn": "MHC_I",
        }
    )
    expected = {"tool": "pvacseq", "method": "cterm", "threshold": 0.5}
    actual = somatic_neoepitope_prediction_workflow.get_args("netchop", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {"threads": 8, "time": "143:59:59", "memory": "32G"}
    # Evaluate
    for resource in expected_dict.keys():
        msg_error = f"Unexpected value of '{resource}' in 'run' sub-step"
        actual = somatic_neoepitope_prediction_workflow.get_resource("netchop", "run", resource)()
        assert actual == expected_dict[resource], msg_error

# ---- netMHCstab

def test_somatic_neoepitope_prediction_netstab_set_part_get_input_files(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
            "tool": "pvacseq",
        }
    )
    annotated_tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "alleles": f"work/{annotated_tpl}/out/{annotated_tpl}.hla_types.txt",
        "epitopes": "work/bwa.mutect2.vep.pvacseq.P001-T1-DNA1-WGS1/out/MHC_Class_I/P001-T1-DNA1-WGS1.MHC_I.netchop.tsv",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("netstab", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_netstab_set_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    expected = {
        "netstab": "work/{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}/out/MHC_Class_I/{tumor_dna}.MHC_I.netstab.tsv",
    }
    actual = somatic_neoepitope_prediction_workflow.get_output_files("netstab", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_netstab_set_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/netstab.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("netstab", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_netstab_step_part_get_args(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
            "tool": "pvacseq",
        }
    )
    expected = {
        "tool": "pvacseq",
        "lengths": [8, 9, 10, 11],
        "exclude_bind": ["container", "alleles"],
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("netstab", "run")(wildcards)
    assert actual == expected

# ---- HLA types

def test_somatic_neoepitope_prediction_hla_types_set_part_get_input_files(
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
    optitype_tpl = "/HLA_TYPING/output/optitype.{library}/out/optitype.{library}.json"
    arcashla_tpl = "/HLA_TYPING/output/star.arcashla.{library}/out/star.arcashla.{library}.json"
    expected = [
        optitype_tpl.format(library="P001-N1-DNA1-WGS1"),
        optitype_tpl.format(library="P001-T1-DNA1-WGS1"),
        optitype_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
        arcashla_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
    ]
    actual = somatic_neoepitope_prediction_workflow.get_input_files("hla_types", "run")(wildcards)
    assert actual == expected

def test_somatic_neoepitope_prediction_hla_types_set_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    expected = {
        "hla_types": "work/{mapper}.{caller}.{annotator}.{tumor_dna}/out/{mapper}.{caller}.{annotator}.{tumor_dna}.hla_types.txt",
    }
    actual = somatic_neoepitope_prediction_workflow.get_output_files("hla_types", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_hla_types_set_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/hla_types.{{tumor_dna}}")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("hla_types", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_hla_types_step_part_get_args(
    somatic_neoepitope_prediction_workflow,
):
    wildcards = Wildcards(
        fromdict={
            "mapper": "bwa",
            "caller": "mutect2",
            "annotator": "vep",
            "tumor_dna": "P001-T1-DNA1-WGS1",
            "tool": "pvacseq",
        }
    )
    expected = {}
    actual = somatic_neoepitope_prediction_workflow.get_args("hla_types", "run")(wildcards)
    assert actual == expected

# ---- Main workflow

def test_somatic_neoepitope_prediction_workflow(somatic_neoepitope_prediction_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["hla_types", "link_out", "netchop", "netstab", "phasing", "proteome", "pvacfuse", "pvacseq", "pvacsplice", "pvactools"]
    actual = list(sorted(somatic_neoepitope_prediction_workflow.sub_steps.keys()))
    assert actual == expected

    base = "{mapper}.{caller}.{annotator}.{{tool}}.{{sample}}".format(
        mapper="bwa", caller="mutect2", annotator="vep"
    )
    samples = [f"P00{i}-T{t}-DNA1-WGS1" for (i, t) in ((1, 1), (2, 1), (2, 2))]
    tools = ("pvacseq", "pvacsplice", "pvacfuse")
    mhc_classes = (("MHC_Class_I", "MHC_I"), ("MHC_Class_II", "MHC_II"))
    log_exts = ("log", "log.md5", "conda_list.txt", "conda_list.txt.md5", "conda_info.txt", "conda_info.txt.md5")
    expected = []

    for tool in tools:
        tool_samples = samples
        # pvacfuse & pvacsplice require RNA sample, absent for P002-T1
        if tool == "pvacfuse" or tool == "pvacsplice":
            tool_samples = (samples[0], samples[2])

        tpl = f"output/{base}/out/combined/{{sample}}.Combined.{{ext}}.tsv"
        expected += [
            tpl.format(sample=sample, tool=tool, ext=ext)
            for sample in tool_samples
            for ext in ("all_epitopes", "all_epitopes.aggregated", "filtered")
        ]

        tpl = f"output/{base}/log/{{tool}}.{{sample}}.{{ext}}"
        expected += [
            tpl.format(sample=sample, tool=tool, ext=ext)
            for sample in tool_samples
           for ext in log_exts
        ]

        # netchop not enabled for pvacfuse
        if tool != "pvacfuse":
            tool_mhc_classes = mhc_classes
            # pvacsplice algorithms are only class II
            if tool == "pvacsplice":
                tool_mhc_classes = (mhc_classes[1],)

            tpl = f"output/{base}/out/{{mhc_d}}/{{sample}}.{{mhc_fn}}.netchop.tsv"
            expected += [
                tpl.format(sample=sample, tool=tool, mhc_d=mhc_d, mhc_fn=mhc_fn)
                for sample in tool_samples
                for (mhc_d, mhc_fn) in tool_mhc_classes
            ]

            tpl = f"output/{base}/log/netchop.{{mhc_d}}_{{mhc_fn}}.{{sample}}.{{ext}}"
            expected += [
                tpl.format(sample=sample, tool=tool, mhc_d=mhc_d, mhc_fn=mhc_fn, ext=ext)
                for sample in tool_samples
                for (mhc_d, mhc_fn) in tool_mhc_classes
                for ext in log_exts
            ]

        # No class I predictions for pvacsplice
        if tool != "pvacfuse" and tool != "pvacsplice":
            tpl = f"output/{base}/out/MHC_Class_I/{{sample}}.MHC_I.netstab.tsv"
            expected += [
                tpl.format(sample=sample, tool=tool)
                for sample in tool_samples
            ]

            tpl = f"output/{base}/log/netstab.{{sample}}.{{ext}}"
            expected += [
                tpl.format(sample=sample, tool=tool, ext=ext)
                for sample in tool_samples
                for ext in log_exts
            ]

    expected = list(sorted(expected))
    actual = list(sorted(somatic_neoepitope_prediction_workflow.get_result_files()))
    assert expected == actual

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
            vep:
                cache_dir: /path/to/dir/cache

          germline_variant_calling:
            path_ngs_mapping: NGS_MAPPING
            tools: [gatk4_hc]
            gatk4_hc:
              num_threads: 8

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
            phasing:
              enabled: true
              path_combine_variants: COMBINE_VARIANTS
            tools_hla_typing:
              dna:
                class_i: optitype
              rna:
                class_i: optitype
                class_ii: arcashla
            pvacseq:
                algorithms: ['MHCflurry','MHCnuggetsI']
                class_ii_epitope_length: [10, 11]
                class_i_epitope_length: []
                net_chop:
                  enabled: true
                  path_netchop: /path/to/netchop.bin
            pvacfuse:
                algorithms: all_class_i
            pvacsplice:
                algorithms: all
                class_ii_epitope_length: [10, 11]
                genes_of_interest_file: /path/to/genes.txt
                run_reference_proteome_similarity: true
                net_chop:
                  enabled: true
                  path_netchop: /path/to/netchop.bin
                peptide_fasta: /path/to/peptides.fa
                
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


def test_somatic_neoepitope_prediction_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    for tool in ("pvacseq", "pvacfuse", "pvacsplice"):
        expected = {
            "filtered": (
                "work/{mapper}.{caller}.{annotator}." + tool + ".{tumor_dna}/"
                + "{mhc_class_d,MHC_Class_II?|combined}/{sample}.{mhc_class_fn,MHC_II?|Combined}.filtered.tsv"
            ),
        }
        actual = somatic_neoepitope_prediction_workflow.get_output_files(tool, tool)
        assert actual == expected


def test_somatic_neoepitope_prediction_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    for tool in ("pvacseq", "pvacfuse", "pvacsplice"):
        expected = get_expected_log_files_dict(
            base_out=(
                "work/{mapper}.{caller}.{annotator}." + tool + ".{tumor_dna}/log/"
                + "{mhc_class_d,MHC_Class_II?|combined}.{sample}.{mhc_class_fn,MHC_II?|Combined}.filtered"
            )
        )
        actual = somatic_neoepitope_prediction_workflow.get_log_file(tool, tool)
        assert actual == expected

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


def test_somatic_neoepitope_prediction_pvacseq_rename_step_part_get_input_files(
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
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "rename")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_combine_step_part_get_input_files(
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
        "annotated": f"work/{annotated_tpl}/out/{annotated_tpl}.renamed.vcf.gz",
        "gene_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.gene.sf",
        "transcript_tpms": f"GENE_EXPRESSION_QUANTIFICATION/output/{expr_tpl}/out/{expr_tpl}.transcript.sf",
    }

    # Get actual
    actual = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "combine")(wildcards)
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
    optitype_tpl = "/HLA_TYPING/output/optitype.{library}/out/optitype.{library}.json"
    arcashla_tpl = "/HLA_TYPING/output/star.arcashla.{library}/out/star.arcashla.{library}.json"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "vcf": f"work/{annotated_tpl}/out/{annotated_tpl}.combined.vcf.gz",
        "phased": f"work/{annotated_tpl}/out/{annotated_tpl}.phased.vcf.gz",
        "alleles": [
            optitype_tpl.format(library="P001-T1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-N1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
            arcashla_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
        ],
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


def test_somatic_neoepitope_prediction_pvacseq_rename_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.renamed.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "rename")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_combine_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = {"vcf": f"work/{tpl}/out/{tpl}.combined.vcf.gz"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("pvacseq", "combine")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {
        "pileup": {"threads": 1, "time": "03:59:59", "memory": "6G"},
        "rename": {"threads": 1, "time": "00:59:59", "memory": "2G"},
        "combine": {"threads": 1, "time": "03:59:59", "memory": "6G"},
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
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/pileup")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "pileup")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_rename_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/rename")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "rename")
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_combine_step_part_get_log_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/combine")
    actual = somatic_neoepitope_prediction_workflow.get_log_file("pvacseq", "combine")
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
        "extra_args": "--adjust-MQ 0 --delta-BQ 30 --max-BQ 60 --max-depth 250 --min-BQ 1 --min-MQ 0",
        "tumor_sample": "P001-T1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "pileup")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_rename_step_part_get_args(
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
        "tumor_sample": "P001-T1",
        "normal_sample": "P001-N1",
        "tumor_library": "P001-T1-DNA1-WGS1",
        "normal_library": "P001-N1-DNA1-WGS1",
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "rename")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacseq_combine_step_part_get_args(
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
        "tumor_sample": "P001-T1",
        "normal_sample": "P001-N1",
        "extra_args": (
            r"--ensembl-id --annotation 'CSQ' --annotation-description-regex "
            r"'^Consequence annotations from Ensembl VEP. Format: (?P<titles>.+)$' "
            r"--annotation-gene-id 'Gene' --annotation-separator '\|' "
            r"--annotation-transcript-id 'Feature' --format 'salmon' "
            "--use-ensembl-version 'none'"
        ),
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "combine")(wildcards)
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
    input = somatic_neoepitope_prediction_workflow.get_input_files("pvacseq", "pvacseq")(wildcards)
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1",
        "normal_sample": "P001-N1",
        "class_i": ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*15:32", "HLA-C*04:03", "HLA-C*04:04"],
        "class_ii": ["DPB1*14:01"],
        "algorithms": "MHCflurry MHCnuggetsI",
        "exclude_bind": ["container", "alleles"],
        "extra_args": (
            "--aggregate-inclusion-binding-threshold 5000 --aggregate-inclusion-count-limit 25 "
            "--anchor-contribution-threshold 0.8 --class-ii-epitope-length '10,11' "
            "--downstream-sequence-length 1000 --expn-val 1.0 --fasta-size 200 "
            "--minimum-fold-change 0.0 --netmhciipan-version '4.1' --normal-cov 25 --normal-vaf 0.02 "
            "--percentile-threshold 1.0 --percentile-threshold-strategy 'conservative' "
            "--tdna-cov 25 --tdna-vaf 0.1 --top-score-metric 'median' --top-score-metric2 'percentile' "
            "--transcript-prioritization-strategy 'mane_select' --trna-cov 2 --trna-vaf 0.25"
        )
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacseq", "pvacseq")(wildcards, input)
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
    optitype_tpl = "/HLA_TYPING/output/optitype.{library}/out/optitype.{library}.json"
    arcashla_tpl = "/HLA_TYPING/output/star.arcashla.{library}/out/star.arcashla.{library}.json"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "fusions": f"SOMATIC_GENE_FUSION_CALLING/output/arriba.P001-T1-RNA1-mRNA_seq1/out/arriba.P001-T1-RNA1-mRNA_seq1.fusions.tsv",
        "alleles": [
            optitype_tpl.format(library="P001-T1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-N1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
            arcashla_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
        ],
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
    input = somatic_neoepitope_prediction_workflow.get_input_files("pvacfuse", "pvacfuse")(wildcards)
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1",
        "class_i": ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*15:32", "HLA-C*04:03", "HLA-C*04:04"],
        "class_ii": ["DPB1*14:01"],
        "algorithms": "all_class_i",
        "exclude_bind": ["container", "alleles"],
        "extra_args": (
            "--aggregate-inclusion-binding-threshold 5000 --aggregate-inclusion-count-limit 25 "
            "--class-i-epitope-length '8,9,10,11' --downstream-sequence-length 1000 --expn-val 1.0 "
            "--fasta-size 200 --netmhciipan-version '4.1' --percentile-threshold 1.0 "
            "--percentile-threshold-strategy 'conservative' --read-support 5 "
            "--top-score-metric 'median' --top-score-metric2 'percentile'"
       )
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacfuse", "pvacfuse")(wildcards, input)
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
        "annotated": f"SOMATIC_VARIANT_ANNOTATION/output/{annotated_tpl}/out/{annotated_tpl}.vcf.gz",
        "bam": f"NGS_MAPPING/output/{rna_tpl}/out/{rna_tpl}.bam",
        "strandedness": f"NGS_MAPPING/output/{rna_tpl}/strandedness/{rna_tpl}.decision.json",
        "reference": "/path/to/ref.fa",
        "features": "/path/to/gencode.gtf",
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
    optitype_tpl = "/HLA_TYPING/output/optitype.{library}/out/optitype.{library}.json"
    arcashla_tpl = "/HLA_TYPING/output/star.arcashla.{library}/out/star.arcashla.{library}.json"
    expected = {
        "container": "work/containers/out/pvactools.sif",
        "junctions": f"work/{annotated_tpl}/out/{annotated_tpl}.junctions.tsv",
        "alleles": [
            optitype_tpl.format(library="P001-T1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-N1-DNA1-WGS1"),
            optitype_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
            arcashla_tpl.format(library="P001-T1-RNA1-mRNA_seq1"),
        ],
        "annotated": f"work/{annotated_tpl}/out/{annotated_tpl}.renamed.vcf.gz",
        "genes": "/path/to/genes.txt",
        "peptides": "/path/to/peptides.fa",
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
    input = somatic_neoepitope_prediction_workflow.get_input_files("pvacsplice", "pvacsplice")(wildcards)
    expected = {
        "n_threads": 1,
        "tumor_sample": "P001-T1",
        "normal_sample": "P001-N1",
        "class_i": ["HLA-A*02:01", "HLA-A*11:01", "HLA-B*15:32", "HLA-C*04:03", "HLA-C*04:04"],
        "class_ii": ["DPB1*14:01"],
        "algorithms": "all",
        "exclude_bind": ["container", "alleles"],
        "extra_args": (
            "--run-reference-proteome-similarity "
            "--aggregate-inclusion-binding-threshold 5000 --aggregate-inclusion-count-limit 25 "
            "--anchor-types 'A,D,NDA' --class-i-epitope-length '8,9,10,11' --class-ii-epitope-length '10,11' "
            "--expn-val 1.0 --fasta-size 200 --junction-score 10 --netmhciipan-version '4.1' "
            "--normal-cov 25 --normal-vaf 0.02 --percentile-threshold 1.0 "
            "--percentile-threshold-strategy 'conservative' --tdna-cov 25 --tdna-vaf 0.1 "
            "--top-score-metric 'median' --top-score-metric2 'percentile' "
            "--transcript-prioritization-strategy 'mane_select' --trna-cov 2 --trna-vaf 0.25 "
            "--variant-distance 100"
        )
    }
    actual = somatic_neoepitope_prediction_workflow.get_args("pvacsplice", "pvacsplice")(wildcards, input)
    assert actual == expected


def test_somatic_neoepitope_prediction_pvacsplice_step_part_get_resource(
    somatic_neoepitope_prediction_workflow,
):
    # Define expected
    expected_dict = {
        "junction": {"threads": 1, "time": "03:59:59", "memory": "6G"},
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
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/phasing")
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
            "mhc_class_d": "MHC_Class_II",
            "sample": "P001-T1",
            "mhc_class_fn": "MHC_II",
        }
    )
    expected = {
        "vcf": "COMBINE_VARIANTS/output/bwa.combined.P001-T1-DNA1-WGS1/out/bwa.combined.P001-T1-DNA1-WGS1.vcf.gz",
        "epitopes": "work/bwa.mutect2.vep.pvacseq.P001-T1-DNA1-WGS1/MHC_Class_II/P001-T1.MHC_II.filtered.tsv",
        "netchop": "/path/to/netchop.bin",
    }
    actual = somatic_neoepitope_prediction_workflow.get_input_files("netchop", "run")(wildcards)
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_output_files(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}"
    expected = {"epitopes": f"work/{tpl}/{{mhc_class_d}}/{{sample}}.{{mhc_class_fn}}.netchop.tsv"}
    actual = somatic_neoepitope_prediction_workflow.get_output_files("netchop", "run")
    assert actual == expected


def test_somatic_neoepitope_prediction_netchop_step_part_get_log_file(
    somatic_neoepitope_prediction_workflow,
):
    tpl = "{mapper}.{caller}.{annotator}.{tool,pvacseq|pvacsplice|pvacfuse}.{tumor_dna}"
    expected = get_expected_log_files_dict(base_out=f"work/{tpl}/log/{{mhc_class_d}}.{{sample}}.{{mhc_class_fn}}.netchop")
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
            "sample": "P001-T1",
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
    expected_dict = {"run": {"threads": 1, "time": "23:59:59", "memory": "32G"}}
    # Evaluate
    for action, resources in expected_dict.items():
        for resource, expected in resources.items():
            msg_error = f"Unexpected value '{expected}' of '{resource}' in '{action}' sub-step"
            actual = somatic_neoepitope_prediction_workflow.get_resource("netchop", action, resource)()
            assert actual == expected, msg_error

# ---- Main workflow

def test_somatic_neoepitope_prediction_workflow(somatic_neoepitope_prediction_workflow):
    """Test simple functionality of the workflow"""
    # Check created sub steps
    expected = ["link_out", "netchop", "phasing", "pvacfuse", "pvacseq", "pvacsplice", "pvactools"]
    actual = list(sorted(somatic_neoepitope_prediction_workflow.sub_steps.keys()))
    assert actual == expected

    log_exts = ("log", "log.md5", "conda_list.txt", "conda_list.txt.md5", "conda_info.txt", "conda_info.txt.md5")

    # Check result file construction
    tpl = "output/{mapper}.{caller}.{annotator}.{tool}.P00{i}-T{t}-DNA1-WGS1/log/{dirname}.P00{i}-T{t}.{filename}.{ext}.{log_ext}"
    expected = [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", tool=tool, dirname=dirname, filename=filename, i=i, t=t, ext=ext, log_ext=log_ext)
        for (tool, dirname, filename, ext) in (
            ("pvacfuse", "MHC_Class_I", "MHC_I", "filtered"),
            ("pvacseq", "MHC_Class_II", "MHC_II", "netchop"),
            ("pvacsplice", "combined", "Combined", "netchop"),
        )
        for (i, t) in ((1, 1), (2, 2))
        for log_ext in log_exts
    ]

    expected += [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", tool="pvacseq", dirname="MHC_Class_II", filename="MHC_II", ext="netchop", i=2, t=1, log_ext=log_ext)
        for log_ext in log_exts
    ]

    tpl = "output/{mapper}.{caller}.{annotator}.{tool}.P00{i}-T{t}-DNA1-WGS1/{dirname}/P00{i}-T{t}.{filename}.{ext}.tsv"
    expected += [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", tool=tool, dirname=dirname, filename=filename, ext=ext, i=i, t=t)
        for (tool, dirname, filename, ext) in (
            ("pvacfuse", "MHC_Class_I", "MHC_I", "filtered"),
            ("pvacseq", "MHC_Class_II", "MHC_II", "netchop"),
            ("pvacsplice", "combined", "Combined", "netchop"),
        )
        for (i, t) in ((1, 1), (2, 2))
    ]

    expected += [
        tpl.format(mapper="bwa", caller="mutect2", annotator="vep", tool="pvacseq", dirname="MHC_Class_II", filename="MHC_II", ext="netchop", i=2, t=1)
    ]

    expected = list(sorted(expected))
    actual = list(sorted(somatic_neoepitope_prediction_workflow.get_result_files()))
    assert expected == actual

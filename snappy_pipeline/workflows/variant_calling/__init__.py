# -*- coding: utf-8 -*-
"""Small variant calling

Calls small variants (SNVs and indels) on the alignments of the ``ngs_mapping`` step.
The ``tool`` field selects the caller.

``gatk4_hc_gvcf``
    GATK 4 HaplotypeCaller in the gVCF workflow: variant discovery per library, ``CombineGVCFs`` per pedigree, then ``GenotypeGVCFs``.
    The output VCF carries physical phasing information.
    The gVCF files are result files, too.

``gatk4_hc_joint``
    GATK 4 HaplotypeCaller with joint calling and direct VCF generation.
    Useful to explore problems with *de novo* calls that the gVCF workflow may introduce.

``gatk3_hc``
    GATK 3 HaplotypeCaller with joint calling.

``gatk3_ug``
    GATK 3 UnifiedGenotyper with joint calling.

``bcftools_call``
    ``bcftools mpileup | bcftools call``.
    Fast, with lower sensitivity and specificity, so suited for baseline calls in benchmarks.

``mutect2``
    Somatic calling with GATK 4 Mutect2 (tumor-normal pair or tumor only, see ``tumor_normal_mode``), with optional contamination estimation and filtration.

``freebayes``
    Candidate calling for ``varlociraptor``: ``freebayes --pooled-continuous`` jointly over the libraries of each entity (the primary DNA libraries of the pedigree with ``group_by: related``, else the library and its matched normal), in parallel over genome windows (``num_jobs``).
    With ``depends_on.kits``, the candidates are restricted to the targets of the libraries' kit, and the windows hold equal numbers of target bases instead of equal lengths, so the jobs take about equally long.
    The calls are candidates, not genotype calls: they carry the tag ``candidates`` and feed ``depends_on.candidates`` of a ``varlociraptor`` task.

``varlociraptor``
    Calling with varlociraptor against a *scenario* of the entity's libraries, on the candidates of ``depends_on.candidates``: a ``freebayes`` task for small variants, or an ``sv_calling_targeted`` or ``sv_calling_wgs`` task for structural variants (the calls are then tagged ``sv``), with the same ``group_by`` and ``library_selection``.
    The candidates of an entity are split into ``num_jobs`` chunks of equal summed depth (``INFO/DP``, so the observation jobs take about equally long), with the records of a breakend event in one chunk; alignment properties per library, observations per library and chunk, calls per chunk, concatenated.
    The scenario is a preset (``scenario``), rendered from the sample sheets for each entity, or a YTE template (``template``):

    ``germline``
        The pedigree (``group_by: related``, the default) or one library: each member with its sex and the ploidy of ``species`` (diploid when the sex is unknown); a member whose parents are both in the task inherits from them, in any constellation (trios, siblings, several generations, uncles and cousins with their parents).
        The events are ``present`` (a member carries the variant) and ``denovo`` (a child carries a variant that its parents do not), which ``present`` excludes; without a member with both parents in the task, ``denovo`` is the impossible event, so one FDR control over ``present`` and ``denovo`` fits every pedigree.
        An unsequenced relative that connects members (a parent of two sequenced siblings, the grandparents of cousins, the other parent of a child with one sequenced parent) joins the scenario as a sample without observations, named after the donor, so inheritance flows through it and the posteriors integrate over its genotype; the VCF holds the observed libraries only.
        A member with a sequenced parent outside the task (one that ``library_selection`` drops) is an error, since varlociraptor models inheritance from both parents only: take both parents into the task or give a ``template``.
        Libraries without a pedigree (``group_by: library``) are called alone.

    ``somatic``
        Each tumor with its matched normal (``tumor_normal_mode``, as for ``mutect2``), with the default scenario of dna-seq-varlociraptor: the events ``somatic_tumor_low``, ``somatic_tumor_medium`` and ``somatic_tumor_high`` (tumor allele frequency below 0.1, from 0.1 to 0.3, from 0.3), ``somatic_normal``, ``germline_het`` and ``germline_hom``.
        The tumor purity is the sample sheet column ``purity`` (a custom field, from 0 to 1), as the contamination of the tumor by the normal; a tumor without one counts as pure, with a warning.
        A tumor without a matched normal keeps the normal in the scenario as an unobserved sample named ``normal``, so the events stay the same; the VCF then has a sample ``normal`` without observations.
        A donor with several tumors in the task (a primary and a relapse) is an error, since the preset calls one tumor with its normal: call each time point in a task of its own (``library_selection``), or give a ``template`` with ``group_by: related``, which then receives the donor's tumors (the entity is the primary tumor) and the normal for a joint scenario such as the tumor/normal/relapse one of the varlociraptor scenario catalog.

    An FFPE library (sample sheet custom field ``isFfpe``, ``Y`` or ``N``) gets varlociraptor's model of the artifacts of an FFPE fragment (``fragment_end_bias``) in both presets.
    The samples of the scenario are the library names, and so are the samples of the VCF, which carries ``GT`` from ``varlociraptor genotype``.
    The calls are all candidates with the posterior probability of each event (``INFO/PROB_<EVENT>``), not a filtered call set: chain a ``variant_filtration`` task with ``tool: varlociraptor`` (FDR control on the events of interest) before any consumer.
    A ``template`` is rendered with ``samples`` (the library table rows of the entity's libraries, as dicts with ``library_name``, ``role``, ``sex``, ``father_name``, ``mother_name``, ``matched_normal_lib`` and the sheet's other columns; for the germline scenario also the unsequenced connecting relatives, with ``observed: false``), ``library_name`` (the entity) and ``settings`` (``species`` and ``purity``); its sample keys must be the library names.

By default (``group_by: related``), the germline callers call each pedigree jointly and produce one result per pedigree, named after the index library.
With ``group_by: library``, they call each library on its own and produce one result per library.
``mutect2`` produces one result per tumor DNA library; ``library_selection`` can only narrow these.
``varlociraptor`` groups like the germline callers with ``scenario: germline`` and like ``mutect2`` with ``scenario: somatic``.
The GATK callers need ``depends_on.dbsnp`` in addition to ``depends_on.reference``.
With ``mutect2``, ``depends_on.panel_of_normals`` optionally names the panel VCF that Mutect2 uses to filter recurrent artefacts: a ``panel_of_normals`` task with ``tool: mutect2``, or an ``external_data`` task with ``{type: models, tags: [pon, mutect2]}`` and the key ``vcf``, for example ``1000g_pon.hg38.vcf.gz`` of the GATK resource bundle.
``depends_on.germline_resource`` optionally names an ``external_data`` task with the population allele frequencies (``{type: variants, tags: [population]}``, key ``vcf``, for example the af-only gnomAD VCF of the GATK resource bundle).
With ``mutect2.contamination.enabled``, ``depends_on.common_variants`` is required: an ``external_data`` task with ``{type: variants, tags: [common]}`` and the key ``vcf``, for example ``small_exac_common_3.vcf`` of the GATK resource bundle.
Windows over the genome (GATK 4 ``PreprocessIntervals``) are processed in parallel; each caller's section has ``window_length`` and ``padding``, the bases added on both sides of a window.
Each caller's section has ``ignore_contigs``, regular expressions for contigs that it skips (by default decoys, HLA, EBV and the hs37d5 extras); a YAML anchor in ``config.yaml`` shares one list between tasks.
The ``FILTER`` IDs of the callers (``LowQual`` of the GATK callers; ``germline``, ``contamination`` and the other filters of ``FilterMutectCalls``) are reserved: no other task may be named like one of them, nor ``PASS``.

Outputs
-------

Below ``tasks/<task>/output/``, for each ``{library_name}``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.
  For ``mutect2`` this is the filtered (PASS) calls.
- ``{library_name}/out/{library_name}.g.vcf.gz`` with the same companion files, for ``gatk4_hc_gvcf`` only.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
The germline callers produce ``variants`` with tags ``germline``, ``snv``, ``indel``.
``mutect2`` produces ``somatic``, ``snv``, ``indel`` for all calls and additionally ``filtered`` for the PASS calls.
``freebayes`` produces ``candidates``, ``snv``, ``indel``.
``varlociraptor`` produces ``germline`` or ``somatic`` (its ``scenario``), ``snv``, ``indel`` and ``varlociraptor``.

Statistics
----------

The enabled statistics are result files below ``{library_name}/report/``, with ``.md5`` files.

``bcftools_stats``
    ``bcftools stats`` per donor of each pedigree: ``bcftools_stats/{library_name}.{donor_library_name}.txt``.
    One job per task also summarizes the statistics of all samples in the result file ``qc_summary/out/qc_summary.tsv`` (one row per sample and variant set: genotype counts, het/hom and Ts/Tv ratios, depth) and renders tables and charts (substitution types, indel length, depth and quality distributions) for the report below ``work/qc_summary/report/html/``.

``bcftools_roh``
    Runs of homozygosity with ``bcftools roh``: ``roh/{library_name}.txt``.
    With ``depends_on.kits``, the regions are the ``targets`` of the kit of the targeted libraries in the VCF, which must share one kit.
    Without it, ROH are called on the whole genome.

``baf_file_generation``
    B-allele fraction per donor as a BigWig file, to confirm structural variants or runs of homozygosity visually.

Examples
--------

Mutect2 with the resources of the GATK resource bundle, one ``external_data`` task per file:

.. code-block:: yaml

    tasks:
      - name: genome
        step: external_data
        config:
          produces: {type: raw, tags: [reference, dna]}
          files: {fasta: /refs/GRCh38/Homo_sapiens_assembly38.fasta}
      - name: bwa_index
        step: external_data
        config:
          produces: {type: index, tags: [bwa, dna]}
          files: {index: /refs/GRCh38/bwa/Homo_sapiens_assembly38}
      - name: mapping
        step: ngs_mapping
        config:
          depends_on: {reads: data_sets, reference: genome, index: bwa_index}
          tool: bwa
          bwa: {}
      - name: gatk_pon
        step: external_data
        config:
          produces: {type: models, tags: [pon, mutect2]}
          files: {vcf: /refs/GRCh38/gatk/1000g_pon.hg38.vcf.gz}
      - name: gnomad
        step: external_data
        config:
          produces: {type: variants, tags: [population]}
          files: {vcf: /refs/GRCh38/gatk/af-only-gnomad.hg38.vcf.gz}
      - name: exac_common
        step: external_data
        config:
          produces: {type: variants, tags: [common]}
          files: {vcf: /refs/GRCh38/gatk/small_exac_common_3.hg38.vcf.gz}
      - name: calling
        step: variant_calling
        config:
          depends_on:
            alignments: mapping
            reference: genome
            panel_of_normals: gatk_pon
            germline_resource: gnomad
            common_variants: exac_common
          tool: mutect2
          mutect2:
            contamination: {enabled: true}

Varlociraptor on tumor/normal pairs of an exome cohort: candidates from freebayes in the kit's
targets, calls against the somatic scenario, then FDR control on the somatic events, whose result
the annotation and the other consumers read:

.. code-block:: yaml

    tasks:
      - name: genome
        step: external_data
        config:
          produces: {type: raw, tags: [reference, dna]}
          files: {fasta: /refs/GRCh38/Homo_sapiens_assembly38.fasta}
      - name: bwa_index
        step: external_data
        config:
          produces: {type: index, tags: [bwa, dna]}
          files: {index: /refs/GRCh38/bwa/Homo_sapiens_assembly38}
      - name: kits
        step: external_data
        config:
          produces: {type: regions, tags: [kits]}
          kits:
            Agilent_SureSelect_V6: {targets: /refs/kits/V6_targets.bed, baits: /refs/kits/V6_baits.bed}
      - name: mapping
        step: ngs_mapping
        config:
          depends_on: {reads: data_sets, reference: genome, index: bwa_index}
          tool: bwa
          bwa: {}
      - name: candidates
        step: variant_calling
        config:
          depends_on: {alignments: mapping, reference: genome, kits: kits}
          tool: freebayes
          freebayes: {}
      - name: calling
        step: variant_calling
        config:
          depends_on: {alignments: mapping, reference: genome, candidates: candidates}
          tool: varlociraptor
          varlociraptor: {scenario: somatic}
      - name: fdr_control
        step: variant_filtration
        config:
          depends_on: {variants: calling}
          tool: varlociraptor
          varlociraptor: {events: [somatic_tumor_high, somatic_tumor_medium]}

For a germline cohort, ``candidates`` gets ``group_by: related``, ``calling`` gets ``scenario: germline`` and ``fdr_control`` the events ``present`` and ``denovo``.
"""

import math
from pathlib import Path
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
    WritePedigreeStepPart,
    _declared_group_by,
)
from snappy_pipeline.workflows.abstract.common import (
    SnakemakeDict,
    SnakemakeDictItemsGenerator,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_calling.model import ScenarioKind, Tool, TumorNormalMode
from snappy_pipeline.workflows.variant_calling.scenario import (
    IDENTIFIER,
    check_scenario,
    germline_scenario,
    is_ffpe,
    render_template,
    somatic_scenario,
    to_yaml,
)

from .model import VariantCalling as VariantCallingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available germline variant callers
VARIANT_CALLERS = (
    "bcftools_call",
    "gatk3_hc",
    "gatk3_ug",
    "gatk4_hc_gvcf",
    "gatk4_hc_joint",
)

#: FILTER IDs in the header of FilterMutectCalls (GATK 4.6.2) and of the GATK germline callers
MUTECT2_FILTER_IDS = frozenset(
    {
        "FAIL",
        "base_qual",
        "clustered_events",
        "contamination",
        "duplicate",
        "fragment",
        "germline",
        "haplotype",
        "low_allele_frac",
        "map_qual",
        "multiallelic",
        "n_ratio",
        "normal_artifact",
        "orientation",
        "panel_of_normals",
        "position",
        "possible_numt",
        "slippage",
        "strand_bias",
        "strict_strand",
        "weak_evidence",
    }
)
GATK_GERMLINE_FILTER_IDS = frozenset({"LowQual"})


class VariantCallingStepPart(BaseStepPart):
    """Base class for germline variant calling step parts

    Variant calling is performed on a per-pedigree level.  The (single) index individual is used
    for naming the output file.
    """

    #: Class available actions
    actions = ("run",)

    result_keys = ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")

    def _get_log_file(self, action: str) -> SnakemakeDict:
        token = "{library_name}"
        return {"log": f"work/{token}/log/{token}.{self.name}_{action}.log"}

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        df = self.parent.build_library_dataframe()
        if df.empty:
            return

        group_by = self.parent.effective_group_by

        # Determine the target libraries based on the wildcard
        if group_by == "related":
            # wildcards.library_name is the index library name; find its related libraries
            df_target = df[
                df["library_name"].isin(self.parent.get_related_libraries(wildcards.library_name))
            ]
        elif group_by == "library":
            # Per-library logic
            df_target = df[df["library_name"] == wildcards.library_name]
        else:
            raise ValueError(f"Unknown group_by: {group_by}")

        if df_target["kind"].eq("cancer").all():
            yield "ped", []
        else:
            if group_by == "related":
                yield (
                    "ped",
                    f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
                )
            else:
                related_group = df_target.iloc[0]["related_group"]
                yield "ped", f"work/write_pedigree.{related_group}/out/{related_group}.ped"

        bams = []
        for _, row in df_target.iterrows():
            infix = row["library_name"]
            aln: ExpectedAlignments = self.parent.get_upstream_paths(
                "alignments", library_name=infix
            )
            bams.append(aln.bam)
        yield "bam", bams

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        token = "{library_name}"
        work_files = {
            "vcf": f"work/{token}/out/{token}.vcf.gz",
            "vcf_md5": f"work/{token}/out/{token}.vcf.gz.md5",
            "vcf_tbi": f"work/{token}/out/{token}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{token}/out/{token}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()


class BcftoolsCallStepPart(VariantCallingStepPart):
    """Germline variant calling with bcftools"""

    #: Step name
    name = "bcftools_call"

    default_resource_usage = ResourceUsage(
        threads=16, runtime="2d", mem=f"{int(3.75 * 1024 * 16)}MB"
    )

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        parent = super()._get_input_files_run(wildcards)
        for k, v in parent.items():
            yield k, v
        yield from self.parent.reference_inputs().items()
        yield "reference_index", self.parent.get_upstream_paths("reference").fai

    def _get_params_run(self, wildcards):
        reference_path = self.parent.get_upstream_paths("reference").fasta
        if "GRCh37" in reference_path or "hg19" in reference_path:
            assembly = "GRCh37"
        elif "GRCh38" in reference_path or "hg38" in reference_path:
            assembly = "GRCh38"
        else:
            assembly = "unknown"

        return {
            "assembly": assembly,
            "ignore_contigs": self.config.bcftools_call.ignore_contigs,
            "window_length": self.config.bcftools_call.window_length,
            "num_threads": self.config.bcftools_call.num_threads,
            "max_depth": self.config.bcftools_call.max_depth,
            "max_indel_depth": self.config.bcftools_call.max_indel_depth,
            "padding": self.config.bcftools_call.padding,
            "annotate": self.config.bcftools_call.annotate,
            "per_sample_mf": self.config.bcftools_call.per_sample_mf,
            "redo_baq": self.config.bcftools_call.redo_baq,
            "variants_only": self.config.bcftools_call.variants_only,
        }


class GatkCallerStepPartBase(VariantCallingStepPart):
    """Base class for GATK v3/v4 variant callers"""

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        """
        Adds the reference & dbsnp input files to the bams & pedigree.
        Note that dpsnp is not needed by combine & genotype sub-steps of GATK4 hc,
        but the previous code was checking for the presence of dbsnp, so
        when was failing when dbsnp is not present anyway.
        """
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        yield from self.parent.vcf_inputs("dbsnp").items()

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config[self.name].num_threads
        mem_per_thread = 5.5
        mem_total = int(mem_per_thread * num_threads + 0.5)
        return ResourceUsage(threads=num_threads, runtime="2d", mem=f"{mem_total}GB")


class Gatk3HaplotypeCallerStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK v3 HaplotypeCaller"""

    #: Step name
    name = "gatk3_hc"

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return {
            "num_threads": self.config.gatk3_hc.num_threads,
            "window_length": self.config.gatk3_hc.window_length,
            "allow_seq_dict_incompatibility": self.config.gatk3_hc.allow_seq_dict_incompatibility,
            "ignore_contigs": self.config.gatk3_hc.ignore_contigs,
            "padding": self.config.gatk3_hc.padding,
            "sample_ploidy": self.config.gatk3_hc.sample_ploidy,
        }


class Gatk3UnifiedGenotyperStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK v3 UnifiedGenotyper"""

    #: Step name
    name = "gatk3_ug"

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return {
            "num_threads": self.config.gatk3_ug.num_threads,
            "window_length": self.config.gatk3_ug.window_length,
            "allow_seq_dict_incompatibility": self.config.gatk3_ug.allow_seq_dict_incompatibility,
            "downsample_to_coverage": self.config.gatk3_ug.downsample_to_coverage,
            "ignore_contigs": self.config.gatk3_ug.ignore_contigs,
            "padding": self.config.gatk3_ug.padding,
            "sample_ploidy": self.config.gatk3_ug.sample_ploidy,
        }


class Gatk4HaplotypeCallerJointStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller doing joint calling per pedigree"""

    name = "gatk4_hc_joint"

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return {
            "window_length": self.config.gatk4_hc_joint.window_length,
            "num_threads": self.config.gatk4_hc_joint.num_threads,
            "allow_seq_dict_incompatibility": self.config.gatk4_hc_joint.allow_seq_dict_incompatibility,
            "ignore_contigs": self.config.gatk4_hc_joint.ignore_contigs,
            "padding": self.config.gatk4_hc_joint.padding,
            "annotation_groups": self.config.gatk4_hc_joint.annotation_groups,
        }


class Gatk4HaplotypeCallerGvcfStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller and gVCF workflow"""

    name = "gatk4_hc_gvcf"

    actions = ("discover", "combine_gvcfs", "genotype")

    extra_arguments_action = "discover"

    intermediate_outputs = {"discover": ("gvcf", "gvcf_tbi")}

    #: The genotype rule links the combined gVCF as ``{library}.g.vcf.gz``, for re-genotyping
    result_keys = {
        "genotype": (
            "vcf",
            "vcf_md5",
            "vcf_tbi",
            "vcf_tbi_md5",
            "gvcf",
            "gvcf_md5",
            "gvcf_tbi",
            "gvcf_tbi_md5",
        )
    }

    @dictify
    def _get_input_files_discover(self, wildcards):
        yield from self.parent.reference_inputs().items()
        yield from self.parent.vcf_inputs("dbsnp").items()
        infix = wildcards.library_name
        alignments: ExpectedAlignments = self.parent.get_upstream_paths(
            "alignments", library_name=infix
        )
        yield "bam", getattr(alignments, "bam", None) or alignments.bam

    @dictify
    def _get_input_files_combine_gvcfs(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        yield from self.parent.reference_inputs().items()

        df = self.parent.build_library_dataframe()
        if df.empty:
            return

        group_by = self.parent.effective_group_by

        # Determine the target libraries based on the wildcard
        if group_by == "related":
            # wildcards.library_name is the index library name; find its related libraries
            df_target = df[
                df["library_name"].isin(self.parent.get_related_libraries(wildcards.library_name))
            ]
        elif group_by == "library":
            # Per-library logic
            df_target = df[df["library_name"] == wildcards.library_name]
        else:
            raise ValueError(f"Unknown group_by: {group_by}")

        if df_target["kind"].eq("cancer").all():
            yield "ped", []
        else:
            if group_by == "related":
                yield (
                    "ped",
                    f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
                )
            else:
                related_group = df_target.iloc[0]["related_group"]
                yield "ped", f"work/write_pedigree.{related_group}/out/{related_group}.ped"

        gvcfs = []
        for _, row in df_target.iterrows():
            infix = f"gatk4_hc_gvcf_discover.{row['library_name']}"
            gvcfs.append(f"work/{infix}/out/{infix}.g.vcf.gz")
        yield "gvcf", gvcfs
        yield "gvcf_tbi", [f"{gvcf}.tbi" for gvcf in gvcfs]

    @dictify
    def _get_input_files_genotype(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield from self.parent.reference_inputs().items()

        infix = f"gatk4_hc_gvcf_combine_gvcfs.{wildcards.library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"

    @dictify
    def _get_output_files_discover(self) -> SnakemakeDictItemsGenerator:
        infix = "gatk4_hc_gvcf_discover.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"

    @dictify
    def _get_output_files_combine_gvcfs(self) -> SnakemakeDictItemsGenerator:
        infix = "gatk4_hc_gvcf_combine_gvcfs.{library_name}"
        yield "gvcf", f"work/{infix}/out/{infix}.g.vcf.gz"
        yield "gvcf_tbi", f"work/{infix}/out/{infix}.g.vcf.gz.tbi"

    @dictify
    def _get_output_files_genotype(self) -> SnakemakeDictItemsGenerator:
        infix = "{library_name}"
        result = {
            "gvcf": f"work/{infix}/out/{infix}.g.vcf.gz",
            "gvcf_md5": f"work/{infix}/out/{infix}.g.vcf.gz.md5",
            "gvcf_tbi": f"work/{infix}/out/{infix}.g.vcf.gz.tbi",
            "gvcf_tbi_md5": f"work/{infix}/out/{infix}.g.vcf.gz.tbi.md5",
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from result.items()

    def _get_params_discover(self, wildcards) -> dict[str, Any]:
        return self._params_gvcf("discover")

    def _get_params_combine_gvcfs(self, wildcards) -> dict[str, Any]:
        return self._params_gvcf("combine_gvcfs")

    def _get_params_genotype(self, wildcards) -> dict[str, Any]:
        return self._params_gvcf("genotype")

    def _params_gvcf(self, action: str) -> dict[str, Any]:
        config = self.config.gatk4_hc_gvcf
        return {
            "window_length": config.window_length,
            "num_threads": config.num_threads,
            "allow_seq_dict_incompatibility": config.allow_seq_dict_incompatibility,
            "ignore_contigs": config.ignore_contigs,
            "padding": config.padding,
            "annotation_groups": getattr(config, f"{action}_annotation_groups"),
        }


class ReportStepPart(BaseStepPart):
    """Base class for the reports on the called variants, one per index library (and donor)"""

    actions = ("run",)

    entity_wildcard = "index_library_name"

    #: Whether we generate per-donor files.
    report_per_donor: bool

    def _get_log_file(self, action: str) -> SnakemakeDict:
        token = "{index_library_name}"
        donor = "{donor_library_name}." if self.report_per_donor else ""
        return {"log": f"work/{token}/log/{token}.{donor}{self.name}_{action}.log"}

    def _result_wildcards(self, entity: str) -> dict[str, list[str]]:
        result = super()._result_wildcards(entity)
        if self.report_per_donor and self.parent.effective_group_by == "related":
            for members in self.parent.related_libraries.values():
                if entity in members:
                    result["donor_library_name"] = members
                    break
        return result


class BcftoolsStatsStepPart(ReportStepPart):
    """VCF statistics computation with ``bcftools stats``.

    Statistics are computed overall and per-sample
    """

    name = "bcftools_stats"

    report_outputs = {"run": {"txt": ReportOutput("bcftools stats", caption="bcftools_stats.rst")}}

    report_per_donor = True

    result_keys = ("txt", "txt_md5")

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="1024MB")

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{index_library_name}/out/{index_library_name}.vcf.gz".format(**wildcards),
        )

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = (
            "work/{index_library_name}/report/bcftools_stats/"
            "{index_library_name}.{donor_library_name}"
        )
        work_files = {key: f"{base_path}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"donor_library_name": wildcards.donor_library_name}


class BcftoolsStatsSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the ``bcftools stats`` of all samples"""

    source = ("bcftools_stats", "run")
    source_keys = {"txt": "txt"}
    format = "bcftools_stats"
    #: One file per sample of the entity's VCF
    label_wildcard = "donor_library_name"

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "bcftools stats", caption="bcftools_stats_summary.rst", htmlindex="index.html"
            )
        }
    }


class BcftoolsRohStepPart(ReportStepPart):
    """ROH calling with `bcftools roh`."""

    name = "bcftools_roh"
    report_per_donor = False
    result_keys = ("txt", "txt_md5")
    default_resource_usage = ResourceUsage(threads=1, runtime="1h", mem="4GB")

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{index_library_name}/out/{index_library_name}.vcf.gz".format(**wildcards),
        )
        yield "path_af_file", self.config.get(self.name).get("path_af_file")
        if self.parent.config.depends_on.kits:
            kit = self.parent.common_kit("kits", self.vcf_libraries(wildcards.index_library_name))
            yield "path_targets", self.parent.get_upstream_paths("kits", library_kit=kit).targets

    def vcf_libraries(self, index_library_name: str) -> list[str]:
        """Return the libraries in the VCF of the index library."""
        if self.parent.effective_group_by == "library":
            return [index_library_name]
        return self.parent.get_related_libraries(index_library_name)

    def _get_params_run(self, wildcards):
        return {
            name: self.config.bcftools_roh.get(name)
            for name in [
                "ignore_homref",
                "skip_indels",
                "rec_rate",
            ]
        }

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = "work/{index_library_name}/report/roh/{index_library_name}"
        for key, ext in ext_names.items():
            yield key, f"{base_path}{ext}"


class BafFileGenerationStepPart(ReportStepPart):
    """Class for computing B allele fraction files.

    One file is generated per sample in the output VCF files.
    """

    #: Step name
    name = "baf_file_generation"

    #: Whether to report results per donor
    report_per_donor = True

    result_keys = ("bw", "bw_md5")

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="1024MB")

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{index_library_name}/out/{index_library_name}.vcf.gz".format(**wildcards),
        )
        yield "reference_index", self.parent.get_upstream_paths("reference").fai

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        base_path = (
            "{index_library_name}/report/baf/{index_library_name}.{donor_library_name,[^\\.]+}.baf"
        )
        ext_names = {"bw": ".bw", "bw_md5": ".bw.md5"}
        work_files = {}
        for key, ext in ext_names.items():
            work_files[key] = f"work/{base_path}{ext}"
        yield from work_files.items()

    def _get_params_run(self, wildcards):
        return {"min_dp": self.config.baf_file_generation.min_dp}


class SomaticVariantCallingStepPart(BaseStepPart):
    @dictify
    def _get_input_files_run(self, wildcards: Wildcards):
        tumor = self.parent.get_upstream_paths("alignments", library_name=wildcards.library_name)

        input_files = {
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
        }

        normal_library = self.get_normal_lib_name(wildcards)
        if normal_library:
            normal = self.parent.get_upstream_paths("alignments", library_name=normal_library)
            input_files.update(
                {
                    "normal_bam": normal.bam,
                    "normal_bai": normal.bai,
                }
            )

        return input_files

    def get_normal_lib_name(self, wildcards):
        """Return the matched normal library of the tumor library, or None"""
        df = self.parent.build_library_dataframe()
        tumor_df = df[df["library_name"] == wildcards.library_name]
        if tumor_df.empty:
            return None
        return tumor_df.iloc[0].get("matched_normal_lib") or None

    def get_tumor_lib_name(self, wildcards):
        return wildcards.library_name


class Mutect2StepPart(SomaticVariantCallingStepPart):
    name = "mutect2"

    actions = (
        "scatter",
        "run",
        "gather",
        "filter",
        "contamination",
        "pileup_normal",
        "pileup_tumor",
    )

    #: The sub-sections ``filtration`` and ``contamination`` hold the arguments of the other
    #: GATK actions
    extra_arguments_action = "run"

    temp_outputs = {"scatter": ("regions",), "run": ("vcf", "vcf_tbi", "stats", "f1r2")}

    intermediate_outputs = {
        "gather": ("vcf", "vcf_tbi", "stats", "orientation"),
        "pileup_normal": ("pileup",),
        "pileup_tumor": ("pileup",),
        "contamination": ("table", "segments"),
    }

    @property
    def resource_usage(self) -> dict[str, ResourceUsage]:
        return {
            "scatter": ResourceUsage(threads=1, runtime="2m", mem="1000MB"),
            "run": ResourceUsage(
                threads=self.config.mutect2.num_cores or 1, runtime="5d", mem="8000MB"
            ),
            "gather": ResourceUsage(threads=1, runtime="4h", mem="32768MB"),
            "filter": ResourceUsage(threads=2, runtime="4h", mem="15872MB"),
            "contamination": ResourceUsage(threads=2, runtime="4h", mem="7680MB"),
            "pileup_normal": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
            "pileup_tumor": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
        }

    def _get_params_scatter(self, wildcards):
        return {
            "ignore_contigs": sorted(set(self.config.mutect2.ignore_contigs)),
            "padding": self.config.mutect2.padding,
        }

    def _get_params_pileup_normal(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "normal_lib_name": self.get_normal_lib_name(wildcards)
        }

    def _get_params_pileup_tumor(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "tumor_lib_name": self.get_tumor_lib_name(wildcards)
        }

    def _get_params_contamination(self, wildcards):
        return {
            "java_options": self.config.mutect2.contamination.java_options,
            "extra_arguments": self.config.mutect2.contamination.extra_arguments,
        }

    def _get_params_run(self, wildcards):
        return {
            "normal_lib_name": self.get_normal_lib_name(wildcards),
            "java_options": self.config.mutect2.java_options,
            "extra_arguments": self.config.mutect2.extra_arguments,
        }

    def _get_params_filter(self, wildcards):
        return self.config.mutect2.filtration.model_dump(by_alias=True)

    def _get_input_files_scatter(self, wildcards):
        return {"fai": self.parent.get_upstream_paths("reference").fai}

    def _get_input_files_run(self, wildcards):
        scatteritem_base_path = (
            "work/{library_name}/out/{library_name}/mutect2par/scatter/{scatteritem}".format(
                **wildcards
            )
        )

        tumor = self.parent.get_upstream_paths("alignments", library_name=wildcards.library_name)
        input_files = {
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
            "region": scatteritem_base_path + ".region.bed",
        }

        tumor_normal_mode = self.config.mutect2.tumor_normal_mode
        if tumor_normal_mode != TumorNormalMode.TUMOR_ONLY:
            normal_library = self.get_normal_lib_name(wildcards)
            if normal_library:
                normal = self.parent.get_upstream_paths("alignments", library_name=normal_library)
                input_files.update(
                    {
                        "normal_bam": normal.bam,
                        "normal_bai": normal.bai,
                    }
                )
            else:
                if tumor_normal_mode == TumorNormalMode.PAIRED:
                    raise ConfigError(
                        f"Normal sample for tumor {wildcards.library_name} required but not found."
                    )

        input_files.update(self.parent.reference_inputs())

        if self.config.depends_on.germline_resource:
            input_files.update(self.parent.vcf_inputs("germline_resource"))
        if self.config.depends_on.panel_of_normals:
            input_files.update(self.parent.vcf_inputs("panel_of_normals"))

        return input_files

    def _get_input_files_gather(self, wildcards):
        gather = self.parent.workflow.globals.get("gather")
        gather = getattr(gather, self.name)
        scatteritem_base_path = (
            "work/{library_name}/out/{library_name}/mutect2par/run/{{scatteritem}}".format(
                **wildcards
            )
        )
        input_files = {
            "vcf": scatteritem_base_path + ".raw.vcf.gz",
            "vcf_tbi": scatteritem_base_path + ".raw.vcf.gz.tbi",
            "stats": scatteritem_base_path + ".raw.vcf.stats",
            "f1r2": scatteritem_base_path + ".raw.f1r2.tar.gz",
        }
        return dict(map(lambda item: (item[0], gather(item[1])), input_files.items()))

    def _get_input_files_filter(self, wildcards):
        base_path = "work/{library_name}/out/{library_name}".format(**wildcards)
        input_files = {
            "raw": base_path + ".raw.vcf.gz",
            "raw_tbi": base_path + ".raw.vcf.gz.tbi",
            "stats": base_path + ".raw.vcf.stats",
            "orientation": base_path + ".raw.read_orientation_model.tar.gz",
            **self.parent.reference_inputs(),
        }
        if self.get_normal_lib_name(wildcards):
            if self.config.mutect2.contamination.enabled:
                input_files["table"] = base_path + ".contamination.tbl"
                input_files["segments"] = base_path + ".segments.tbl"
        return input_files

    def _get_input_files_pileup_normal(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=self.get_normal_lib_name(wildcards)
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            **self.parent.vcf_inputs("common_variants"),
        }

    def _get_input_files_pileup_tumor(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            **self.parent.vcf_inputs("common_variants"),
        }

    def _get_input_files_contamination(self, wildcards: Wildcards):
        base_path = "work/{library_name}/out/{library_name}".format(**wildcards)
        return {
            "normal": base_path + ".normal.pileup",
            "tumor": base_path + ".tumor.pileup",
            **self.parent.reference_inputs(),
        }

    def _get_output_files_scatter(self):
        scatter = getattr(self.parent.workflow.globals.get("scatter"), self.name)
        template = (
            "work/{{library_name}}/out/{{library_name}}/mutect2par/scatter/{scatteritem}.region.bed"
        )
        return {"regions": scatter(template)}

    def _get_output_files_run(self):
        prefix = "work/{library_name}/out/{library_name}/mutect2par/run/{scatteritem}"
        return {
            "vcf": prefix + ".raw.vcf.gz",
            "vcf_tbi": prefix + ".raw.vcf.gz.tbi",
            "stats": prefix + ".raw.vcf.stats",
            "f1r2": prefix + ".raw.f1r2.tar.gz",
        }

    def _get_output_files_gather(self):
        return self._with_exts(
            {
                "vcf": ".raw.vcf.gz",
                "vcf_tbi": ".raw.vcf.gz.tbi",
                "stats": ".raw.vcf.stats",
                "orientation": ".raw.read_orientation_model.tar.gz",
            }
        )

    def _get_output_files_pileup_normal(self):
        return self._with_exts({"pileup": ".normal.pileup"})

    def _get_output_files_pileup_tumor(self):
        return self._with_exts({"pileup": ".tumor.pileup"})

    def _get_output_files_contamination(self):
        return self._with_exts({"table": ".contamination.tbl", "segments": ".segments.tbl"})

    def _get_output_files_filter(self):
        return self._with_exts(
            {
                "vcf": ".vcf.gz",
                "vcf_md5": ".vcf.gz.md5",
                "vcf_tbi": ".vcf.gz.tbi",
                "vcf_tbi_md5": ".vcf.gz.tbi.md5",
                "full_vcf": ".full.vcf.gz",
                "full_vcf_tbi": ".full.vcf.gz.tbi",
            }
        )

    @staticmethod
    def _with_exts(exts: dict[str, str]) -> dict[str, str]:
        return {key: "work/{library_name}/out/{library_name}" + ext for key, ext in exts.items()}

    def _get_log_file(self, action):
        if action == "run":
            postfix = ".{scatteritem}"
        elif action == "gather":
            postfix = ""
        else:
            postfix = "." + action
        return {"log": f"work/{{library_name}}/log/{{library_name}}{postfix}.log"}

    @listify
    def get_result_files(self):
        tpl = "output/{library_name}/out/{library_name}.vcf.gz"
        for library_name in self.parent.libraries["library_name"]:
            for ext in ("", ".md5", ".tbi", ".tbi.md5"):
                yield tpl.format(library_name=library_name) + ext


class ScatteredCallerStepPart(BaseStepPart):
    """Base class of the callers that call an entity in scattered parts and concatenate them
    with the ``bcftools/concat`` wrapper."""

    #: The calls of one part of an entity, a BCF with its ``.csi``
    chunk_tpl = "work/{library_name}/out/{library_name}/calls/{scatteritem}"

    def _scattergather(self, kind: str):
        return getattr(self.parent.workflow.globals.get(kind), self.name)

    def members(self, library_name: str) -> list[str]:
        """Return the libraries called together for the entity ``library_name``, the entity
        first: with ``group_by: related`` the primary DNA libraries of its related group (one
        per individual), else the library and its matched normal, if it has one."""
        if self.parent.effective_group_by == "related":
            primary = set(self.parent.primary_dna_libraries()["library_name"])
            others = [
                name
                for name in self.parent.get_related_libraries(library_name)
                if name in primary and name != library_name
            ]
            return [library_name, *others]
        normal = self.parent.library_row(library_name).get("matched_normal_lib") or None
        return [library_name, normal] if normal else [library_name]

    def _get_input_files_gather(self, wildcards):
        gather = self._scattergather("gather")
        chunk = self.chunk_tpl.replace("{scatteritem}", "{{scatteritem}}").format(**wildcards)
        return {"vcf": gather(chunk + ".bcf"), "vcf_csi": gather(chunk + ".bcf.csi")}

    def _get_output_files_gather(self):
        token = "{library_name}"
        return {
            "vcf": f"work/{token}/out/{token}.vcf.gz",
            "vcf_md5": f"work/{token}/out/{token}.vcf.gz.md5",
            "vcf_tbi": f"work/{token}/out/{token}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{token}/out/{token}.vcf.gz.tbi.md5",
        }

    def _get_log_file_gather(self):
        return {"log": "work/{library_name}/log/{library_name}.gather.log"}


class FreebayesStepPart(ScatteredCallerStepPart):
    """Candidate calling with freebayes, jointly over the libraries of an entity, per genome
    window; the windows are shared by the task's entities"""

    name = "freebayes"

    #: The windows of the task; ``padding`` 0 keeps them disjoint, so no call is made twice
    region_tpl = "work/scatter/{scatteritem}.region.bed"

    actions = ("scatter", "run", "gather")

    extra_arguments_action = "run"

    temp_outputs = {"scatter": ("regions",), "run": ("bcf", "bcf_csi")}

    result_keys = {"gather": ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")}

    resource_usage = {
        "scatter": ResourceUsage(threads=1, runtime="2m", mem="1000MB"),
        "run": ResourceUsage(threads=1, runtime="1d", mem="8GB"),
        "gather": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
    }

    def _get_input_files_scatter(self, wildcards):
        """The genome index and, with ``depends_on.kits``, the targets of the kits of the
        task's libraries: the windows then hold equal target bases instead of equal lengths"""
        inputs = {"fai": self.parent.get_upstream_paths("reference").fai}
        if self.config.depends_on.kits:
            kits = {
                self.parent.library_kits[name]
                for entity in self.parent.output_entities
                for name in self.members(entity)
                if name in self.parent.library_kits
            }
            inputs["targets"] = [
                self.parent.get_upstream_paths("kits", library_kit=kit).targets
                for kit in sorted(kits)
            ]
        return inputs

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from self.parent.reference_inputs().items()
        members = self.members(wildcards.library_name)
        alignments = [
            self.parent.get_upstream_paths("alignments", library_name=name) for name in members
        ]
        yield "bam", [aln.bam for aln in alignments]
        yield "bai", [aln.bai for aln in alignments]
        yield "region", self.region_tpl.format(**wildcards)
        if self.config.depends_on.kits:
            kit = self.parent.common_kit("kits", members)
            yield "targets", self.parent.get_upstream_paths("kits", library_kit=kit).targets

    def _get_params_scatter(self, wildcards):
        return {"ignore_contigs": sorted(set(self.config.freebayes.ignore_contigs)), "padding": 0}

    def _get_params_run(self, wildcards):
        return {
            "min_alternate_fraction": self.config.freebayes.min_alternate_fraction,
            "min_alternate_count": self.config.freebayes.min_alternate_count,
        }

    def _get_output_files_scatter(self):
        return {"regions": self._scattergather("scatter")(self.region_tpl)}

    def _get_output_files_run(self):
        return {"bcf": self.chunk_tpl + ".bcf", "bcf_csi": self.chunk_tpl + ".bcf.csi"}

    def _get_log_file_scatter(self):
        return {"log": "work/scatter/log.log"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.run.{scatteritem}.log"}


def _connecting_donors(rows: dict[str, dict[str, Any]], observed: set[str]) -> set[str]:
    """Return the unsequenced donors of a pedigree that connect observed donors: a donor with an
    observed descendant that has an observed or connecting ancestor, that shares a child with
    an observed or connecting co-parent, or whose children lead to observed donors on two lines.
    The result grows until it is stable."""
    children: dict[str, list[str]] = {}
    for name, row in rows.items():
        for parent in (row["father_name"], row["mother_name"]):
            if parent in rows:
                children.setdefault(parent, []).append(name)

    def leads_to_observed(name: str) -> bool:
        return name in observed or any(leads_to_observed(child) for child in children.get(name, []))

    def ancestors(name: str):
        for parent in (rows[name]["father_name"], rows[name]["mother_name"]):
            if parent in rows:
                yield parent
                yield from ancestors(parent)

    included: set[str] = set()
    candidates = [
        name
        for name, row in rows.items()
        if not row["library_name"] and name not in observed and leads_to_observed(name)
    ]
    while True:
        connected = observed | included
        added = False
        for name in candidates:
            if name in included:
                continue
            co_parents = {
                parent
                for child in children.get(name, [])
                for parent in (rows[child]["father_name"], rows[child]["mother_name"])
                if parent != name and parent in rows
            }
            lines = sum(leads_to_observed(child) for child in children.get(name, []))
            if (
                any(ancestor in connected for ancestor in ancestors(name))
                or co_parents & connected
                or lines >= 2
            ):
                included.add(name)
                added = True
        if not added:
            return included


def _builtin(value: Any) -> Any:
    """Return a library table value as a builtin: numpy scalars unwrapped, NaN as ``None``."""
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, float) and math.isnan(value):
        return None
    return value


class VarlociraptorStepPart(ScatteredCallerStepPart):
    """Calling with varlociraptor on the candidates of ``depends_on.candidates``: the candidates
    of an entity split into chunks of equal depth, alignment properties per library,
    observations per library and chunk, and calls per chunk against the scenario of the entity"""

    name = "varlociraptor"

    actions = ("scatter", "alignment_properties", "preprocess", "call", "gather")

    #: The ``preprocess`` section passes its own extra arguments
    extra_arguments_action = "call"

    temp_outputs = {
        "scatter": ("candidates",),
        "preprocess": ("bcf",),
        "call": ("bcf", "bcf_csi"),
    }

    intermediate_outputs = {"alignment_properties": ("json",)}

    result_keys = {"gather": ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")}

    resource_usage = {
        "scatter": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
        "alignment_properties": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
        "preprocess": ResourceUsage(threads=1, runtime="1d", mem="4GB"),
        "call": ResourceUsage(threads=1, runtime="1d", mem="4GB"),
        "gather": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
    }

    #: The alignment properties of a library, shared by the entities that observe it
    properties_tpl = "work/{member}/out/{member}.alignment_properties.json"

    #: One chunk of the candidates of an entity
    candidates_tpl = "work/{library_name}/out/{library_name}/candidates/{scatteritem}.bcf"

    #: The observations of a library at one chunk of the candidates of an entity
    observations_tpl = (
        "work/{library_name}/out/{library_name}/observations/{member}.{scatteritem}.bcf"
    )

    def __init__(self, parent):
        super().__init__(parent)
        #: The rendered scenario of each entity, built once
        self._scenarios: dict[str, str] = {}

    def members(self, library_name: str) -> list[str]:
        """The observed libraries of the entity: the pedigree or the library alone for the
        germline scenario; for the somatic one the tumor, with ``group_by: related`` the other
        tumors of its donor that the task takes, and its matched normal by
        ``tumor_normal_mode``."""
        cfg = self.config.varlociraptor
        if cfg.scenario == ScenarioKind.germline:
            if self.parent.effective_group_by == "related":
                return super().members(library_name)
            return [library_name]
        members = [library_name]
        if self.parent.effective_group_by == "related":
            members += self.other_tumors(library_name)
        if cfg.tumor_normal_mode != TumorNormalMode.TUMOR_ONLY:
            normal = self.parent.library_row(library_name).get("matched_normal_lib") or None
            if normal:
                members.append(normal)
        return members

    def other_tumors(self, library_name: str) -> list[str]:
        """The other tumor libraries of the donor of ``library_name`` that the task takes"""
        libraries = self.parent.libraries
        donor = self.parent.library_row(library_name)["related_group"]
        same_donor = libraries["related_group"] == donor
        other = libraries["library_name"] != library_name
        return libraries[same_donor & other]["library_name"].tolist()

    def sample(self, library_name: str) -> dict[str, Any]:
        """The library table row of a library with builtin values, as a scenario template sees it"""
        return {
            key: _builtin(value) for key, value in self.parent.library_row(library_name).items()
        }

    def pedigree_samples(self, library_name: str) -> list[dict[str, Any]]:
        """The samples of the germline scenario of a pedigree: its observed members, and the
        unsequenced donors that connect members (a parent of two sequenced children, the
        grandparents of cousins, the other parent of a child with one sequenced parent) as
        samples without observations, named after the donor. Parents are the sample names, or
        ``"0"``; a sequenced parent outside the task keeps its library name, so that the preset
        rejects it."""
        members = self.members(library_name)
        samples = [self.sample(name) for name in members]
        if self.parent.effective_group_by != "related":
            return samples
        donors = self.parent.pedigree_donors(library_name)
        if donors.empty:
            return samples
        rows = {row["donor_name"]: row for row in donors.to_dict("records")}
        by_library = {
            row["library_name"]: name for name, row in rows.items() if row["library_name"]
        }
        observed = {by_library[name] for name in members if name in by_library}
        included = _connecting_donors(rows, observed)

        def parent_name(donor: str) -> str:
            if donor in observed or donor not in rows:
                return rows[donor]["library_name"] if donor in rows else "0"
            if donor in included:
                return donor
            return rows[donor]["library_name"] or "0"

        for sample in samples:
            donor = by_library.get(sample["library_name"])
            if donor is not None:
                sample["father_name"] = parent_name(rows[donor]["father_name"])
                sample["mother_name"] = parent_name(rows[donor]["mother_name"])
        for donor in sorted(included):
            samples.append(
                {
                    "library_name": donor,
                    "sex": rows[donor]["sex"],
                    "father_name": parent_name(rows[donor]["father_name"]),
                    "mother_name": parent_name(rows[donor]["mother_name"]),
                    "observed": False,
                }
            )
        return samples

    def purity(self, library_name: str) -> float | None:
        """The tumor purity of a library from the sample sheet column ``purity``, or ``None``"""
        value = self.sample(library_name).get("purity")
        if value is None or value == "":
            return None
        try:
            purity = float(value)
        except ValueError:
            purity = math.nan
        if not 0 <= purity <= 1:
            raise ConfigError(
                f"the purity of {library_name} is {value!r}, but a number from 0 to 1 is needed"
            )
        return purity

    def scenario(self, library_name: str) -> str:
        """The scenario YAML of the entity, rendered once per run"""
        if library_name not in self._scenarios:
            self._scenarios[library_name] = to_yaml(self._build_scenario(library_name))
        return self._scenarios[library_name]

    def _build_scenario(self, library_name: str) -> dict[str, Any]:
        cfg = self.config.varlociraptor
        members = self.members(library_name)
        somatic = cfg.scenario == ScenarioKind.somatic
        purity = self.purity(library_name) if somatic else None
        if cfg.template:
            settings = {"species": cfg.species.model_dump(), "purity": purity}
            samples = (
                self.pedigree_samples(library_name)
                if not somatic
                else [self.sample(name) for name in members]
            )
            scenario = render_template(
                Path(cfg.template).read_text(), samples, library_name, settings
            )
        elif somatic:
            if others := self.other_tumors(library_name):
                raise ConfigError(
                    "the somatic scenario calls one tumor with its matched normal, but the donor "
                    f"of {library_name} has the tumors {', '.join([library_name, *others])} in "
                    "the task; call each time point in a task of its own (library_selection), or "
                    "give a scenario template with group_by: related for a joint scenario of the "
                    "donor's tumors"
                )
            normal = members[1] if len(members) > 1 else None
            ffpe = frozenset(name for name in members if is_ffpe(self.sample(name).get("is_ffpe")))
            scenario = somatic_scenario(library_name, normal, purity, ffpe)
        else:
            scenario = germline_scenario(
                self.pedigree_samples(library_name), cfg.species.model_dump()
            )
        check_scenario(scenario, members)
        return scenario

    def _get_input_files_scatter(self, wildcards):
        candidates = self.parent.get_upstream_paths(
            "candidates", library_name=wildcards.library_name
        )
        return {"candidates": candidates.vcf}

    def _get_output_files_scatter(self):
        template = self.candidates_tpl.replace("{library_name}", "{{library_name}}")
        return {"candidates": self._scattergather("scatter")(template)}

    def _get_input_files_alignment_properties(self, wildcards):
        alignments = self.parent.get_upstream_paths("alignments", library_name=wildcards.member)
        return {**self.parent.reference_inputs(), "bam": alignments.bam, "bai": alignments.bai}

    def _get_output_files_alignment_properties(self):
        return {"json": self.properties_tpl}

    def _get_input_files_preprocess(self, wildcards):
        alignments = self.parent.get_upstream_paths("alignments", library_name=wildcards.member)
        return {
            **self.parent.reference_inputs(),
            "candidates": self.candidates_tpl.format(**wildcards),
            "alignment_properties": self.properties_tpl.format(**wildcards),
            "bam": alignments.bam,
            "bai": alignments.bai,
        }

    def _get_params_preprocess(self, wildcards):
        preprocess = self.config.varlociraptor.preprocess
        return {"max_depth": preprocess.max_depth, "extra_arguments": preprocess.extra_arguments}

    def _get_output_files_preprocess(self):
        return {"bcf": self.observations_tpl}

    def _get_input_files_call(self, wildcards):
        return {
            "observations": [
                self.observations_tpl.format(
                    library_name=wildcards.library_name,
                    member=member,
                    scatteritem=wildcards.scatteritem,
                )
                for member in self.members(wildcards.library_name)
            ]
        }

    def _get_params_call(self, wildcards):
        """The sample names of the observations, in their order, and the scenario"""
        return {
            "obs_names": self.members(wildcards.library_name),
            "scenario": self.scenario(wildcards.library_name),
        }

    def _get_output_files_call(self):
        return {"bcf": self.chunk_tpl + ".bcf", "bcf_csi": self.chunk_tpl + ".bcf.csi"}

    def _get_log_file_scatter(self):
        return {"log": "work/{library_name}/log/{library_name}.scatter.log"}

    def _get_log_file_alignment_properties(self):
        return {"log": "work/{member}/log/{member}.alignment_properties.log"}

    def _get_log_file_preprocess(self):
        return {
            "log": "work/{library_name}/log/{library_name}.preprocess.{member}.{scatteritem}.log"
        }

    def _get_log_file_call(self):
        return {"log": "work/{library_name}/log/{library_name}.call.{scatteritem}.log"}


class VariantCallingWorkflow(BaseStep):
    """Workflow implementation for germline variant calling"""

    name = "variant_calling"
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "filtered"})),
        DataSignature(DataType.VARIANTS, frozenset({"candidates", "snv", "indel"})),
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel", "varlociraptor"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "varlociraptor"})),
        DataSignature(DataType.VARIANTS, frozenset({"germline", "sv", "varlociraptor"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "sv", "varlociraptor"})),
    ]
    config_model_class = VariantCallingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """Germline small variants; for mutect2 all somatic calls plus the PASS calls; for
        freebayes candidates; for varlociraptor the calls of its scenario, tagged with the tool
        for the FDR control of ``variant_filtration``."""
        if config.tool == Tool.mutect2:
            calls = DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"}))
            return (calls, calls.with_tags("filtered"))
        if config.tool == Tool.freebayes:
            return (DataSignature(DataType.VARIANTS, frozenset({"candidates", "snv", "indel"})),)
        if config.tool == Tool.varlociraptor:
            # The calls are of the kinds of the candidates: structural variants, or small ones
            kinds = {"sv"} if "sv" in upstream["candidates"].tags else {"snv", "indel"}
            tags = {str(config.varlociraptor.scenario), "varlociraptor", *kinds}
            return (DataSignature(DataType.VARIANTS, frozenset(tags)),)
        return (DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel"})),)

    @classmethod
    def default_group_by(cls, config):
        """The germline callers and the germline scenario of varlociraptor call related
        libraries (each pedigree) jointly; mutect2, freebayes and the somatic scenario call each
        library (with its matched normal)."""
        if config.tool == Tool.varlociraptor:
            return "related" if config.varlociraptor.scenario == ScenarioKind.germline else None
        return "related" if config.tool in VARIANT_CALLERS else None

    @classmethod
    def filter_ids(cls, config):
        """The FILTER IDs of mutect2 and of the GATK callers; the other callers write none."""
        if config.tool == Tool.mutect2:
            return MUTECT2_FILTER_IDS
        if config.tool in (Tool.bcftools_call, Tool.freebayes, Tool.varlociraptor):
            return frozenset()
        return GATK_GERMLINE_FILTER_IDS

    #: Default relationship for somatic (mutect2) callers: resolve matched normal.
    default_relationships = {
        "matched_normal_lib": RelationshipDefinition(
            via="donor_name",
            target="role == 'normal' and extraction_type == 'dna'",
        )
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local VCF output paths for a germline-variants signature.

        Arguments:
            signature: Expected to satisfy
                ``DataSignature(DataType.VARIANTS, frozenset({"germline"}))``.
            **kwargs: Accepts ``library_name`` for concrete path rendering; falls back to
                the ``{library_name}`` wildcard placeholder.
        """
        lib = kwargs.get("library_name", "{library_name}")
        # mutect2 writes all calls to .full.vcf.gz and the PASS calls to .vcf.gz
        unfiltered = (
            config.tool == Tool.mutect2
            and signature is not None
            and "filtered" not in signature.tags
        )
        prefix = f"work/{lib}/out/{lib}.full" if unfiltered else f"work/{lib}/out/{lib}"
        return {"vcf": f"{prefix}.vcf.gz", "vcf_tbi": f"{prefix}.vcf.gz.tbi"}

    @classmethod
    def task_library_requirement(cls, config):
        """mutect2 and the somatic scenario of varlociraptor call tumor DNA libraries, which
        need a matched normal or none; ``tumor_normal_mode: paired`` leaves out the tumors
        without one."""
        tumor = "role == 'tumor' and extraction_type == 'dna'"
        if config.tool == Tool.mutect2:
            return tumor
        if config.tool == Tool.varlociraptor:
            cfg = config.varlociraptor
            if cfg.scenario != ScenarioKind.somatic:
                return None
            if cfg.tumor_normal_mode == TumorNormalMode.PAIRED:
                return MATCHED_NORMAL_REQUIREMENT
            return tumor
        return None

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        match self.config.tool:
            case "bcftools_call":
                selected_caller = BcftoolsCallStepPart
            case "gatk3_hc":
                selected_caller = Gatk3HaplotypeCallerStepPart
            case "gatk3_ug":
                selected_caller = Gatk3UnifiedGenotyperStepPart
            case "gatk4_hc_joint":
                selected_caller = Gatk4HaplotypeCallerJointStepPart
            case "gatk4_hc_gvcf":
                selected_caller = Gatk4HaplotypeCallerGvcfStepPart
            case "mutect2":
                selected_caller = Mutect2StepPart
            case "freebayes":
                selected_caller = FreebayesStepPart
            case "varlociraptor":
                selected_caller = VarlociraptorStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes(
            (
                WritePedigreeStepPart,
                selected_caller,
                BcftoolsStatsStepPart,
                BcftoolsRohStepPart,
                BafFileGenerationStepPart,
            )
        )
        if self.config.bcftools_stats and self.config.bcftools_stats.enabled:
            self.register_sub_step_classes((BcftoolsStatsSummaryStepPart,))
        if (
            self.config.depends_on.kits
            and self.config.bcftools_roh
            and self.config.bcftools_roh.enabled
        ):
            for index_library_name in self.output_entities:
                self.common_kit(
                    "kits", self.sub_steps["bcftools_roh"].vcf_libraries(index_library_name)
                )
        if self.config.tool == Tool.freebayes and self.config.depends_on.kits:
            # The candidates of an entity are restricted to the targets of its libraries' kit
            for library_name in self.output_entities:
                self.common_kit("kits", self.sub_steps["freebayes"].members(library_name))
        if self.config.tool == Tool.varlociraptor:
            self._check_varlociraptor_setup()

    def _check_varlociraptor_setup(self) -> None:
        """The candidates task writes its calls per entity of this task, and the library names
        are sample names that a scenario accepts."""
        candidates = self.resolve_dependency("candidates").task_name
        upstream_group_by = _declared_group_by(self.project, candidates, frozenset()) or "library"
        if upstream_group_by != self.effective_group_by:
            raise ConfigError(
                f"the task calls with group_by: {self.effective_group_by}, but its candidates "
                f"task {candidates!r} writes them with group_by: {upstream_group_by}"
            )
        part = self.sub_steps["varlociraptor"]
        invalid = sorted(
            name
            for entity in self.output_entities
            for name in part.members(entity)
            if not IDENTIFIER.match(name)
        )
        if invalid:
            raise ConfigError(
                "varlociraptor needs sample names of letters, digits, '_', '.' and '-', unlike "
                + ", ".join(invalid)
            )
        # Render every scenario now, so that a constellation the preset or the template does
        # not model fails before Snakemake builds the jobs
        for entity in self.output_entities:
            part.scenario(entity)

    def library_warnings(self) -> list[str]:
        """The tumors without a matched normal that ``paired`` mode leaves out, and the tumors
        that the somatic scenario of varlociraptor takes as pure for lack of a purity."""
        warnings = super().library_warnings()
        if (
            self.config.tool == Tool.varlociraptor
            and self.config.varlociraptor.scenario == ScenarioKind.somatic
        ):
            part = self.sub_steps["varlociraptor"]
            unknown = [name for name in self.output_entities if part.purity(name) is None]
            if unknown:
                warnings.append(
                    f"Task {self.task_name!r} takes the tumors without a purity as pure, with no "
                    f"contamination by the normal: {', '.join(unknown)}; the sample sheet column "
                    "purity (a custom field, from 0 to 1) sets it, from the lab, the collaborator "
                    "or a purity estimate (somatic_purity_ploidy_estimate)"
                )
        return warnings

    def get_result_files(self) -> list[str]:
        """The files of the caller and of the enabled reports"""
        names = [self.config.tool]
        for name in ("baf_file_generation", "bcftools_stats", "bcftools_roh"):
            name_config = self.config.get(name)
            if name_config and name_config.enabled:
                names.append(name)
        if "qc_summary" in self.sub_steps:  # the summary of the bcftools stats
            names.append("qc_summary")
        return [path for name in names for path in self.sub_steps[name].get_result_files()]

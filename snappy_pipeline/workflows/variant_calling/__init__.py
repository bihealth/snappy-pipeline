# -*- coding: utf-8 -*-
"""Small variant calling

The step calls small variants (SNVs and indels) on the alignments of the ``ngs_mapping`` step.
The ``tool`` field selects the caller.

``gatk4_hc_gvcf``
    GATK 4 HaplotypeCaller in the gVCF workflow: variant discovery per library, ``CombineGVCFs`` per pedigree, then ``GenotypeGVCFs``.
    The output VCF carries physical phasing information.
    The gVCF files are result files, too.

``gatk4_hc_joint``
    GATK 4 HaplotypeCaller with joint calling and direct VCF generation.
    It helps to explore problems with *de novo* calls that the gVCF workflow may introduce.

``gatk3_hc``
    GATK 3 HaplotypeCaller with joint calling.

``gatk3_ug``
    GATK 3 UnifiedGenotyper with joint calling.

``bcftools_call``
    ``bcftools mpileup | bcftools call``.
    It is fast but has lower sensitivity and specificity, which suits baseline calls in benchmarks.
    ``assembly`` (``GRCh37`` or ``GRCh38``) gives the ploidy of the sex chromosomes.

``mutect2``
    Somatic calling with GATK 4 Mutect2 (tumor-normal pair or tumor only, see ``tumor_normal_mode``), with optional contamination estimation and filtration.
    The genome is split once per task into ``num_jobs`` regions, padded by ``padding``, and each tumor is called in each region.

``freebayes``
    Candidate calling for ``varlociraptor``: ``freebayes --pooled-continuous`` jointly over the libraries of each entity, in parallel over genome windows (``num_jobs``).
    These libraries are the primary DNA libraries of the pedigree with ``group_by: related``, else the library and its matched normal.
    With ``depends_on.kits``, the candidates are restricted to the targets of the libraries' kit.
    The windows then hold equal numbers of target bases instead of equal lengths, so the jobs take about equally long.
    The candidates are normalised against the reference: left-aligned indels, one record per allele, exact duplicates removed.
    With ``atomize`` (the default), complex alleles are split into their SNVs and indels, as the dna-seq-varlociraptor workflow does.
    The calls are candidates, not genotype calls.
    They carry the tag ``candidates`` and feed ``depends_on.candidates`` of a ``varlociraptor`` task.

``varlociraptor``
    Calling with varlociraptor against a *scenario* of the entity's libraries, on the candidates of ``depends_on.candidates``.
    The candidates come from a ``freebayes`` task for small variants, or from an ``sv_calling_targeted`` or ``sv_calling_wgs`` task for structural variants.
    With structural variants, the calls are tagged ``sv``.
    The candidates task must have the same ``group_by`` and ``library_selection`` as this task.
    The candidates of an entity are split into ``num_jobs`` chunks of equal summed depth (``INFO/DP``), so the observation jobs take about equally long.
    The records of a breakend event stay in one chunk.
    The task computes the alignment properties per library, the observations per library and chunk, and the calls per chunk, and concatenates the calls.
    The scenario is a preset (``scenario``), rendered from the sample sheets for each entity, or a YTE template (``template``):

    ``germline``
        The scenario holds the pedigree (``group_by: related``, the default) or one library.
        Each member has its sex and the ploidy of ``species``, and a member of unknown sex is diploid.
        A member whose parents are both in the task inherits from them.
        This works in any constellation: trios, siblings, several generations, uncles and cousins with their parents.
        An unsequenced relative that connects members joins the scenario as a sample without observations, named after the donor.
        Examples are a parent of two sequenced siblings and the grandparents of cousins.
        Inheritance flows through such a relative, and the posteriors integrate over its genotype.
        Varlociraptor models inheritance from both parents only.
        A member with one parent in the scenario and an unknown other parent therefore gets an unobserved second parent, named ``<donor>.father`` or ``<donor>.mother``.
        The VCF holds the observed libraries only.
        The events partition the configurations in which a member carries the variant, by the sheet's ``isAffected``.
        ``denovo_<member>`` exists for each child whose parents are both in the scenario: the child carries a variant that its parents do not.
        Such a variant is a de novo mutation in that member, which its descendants inherit.
        ``<member>`` is the donor name in lower case.
        ``recessive`` means that every affected member is homozygous and no unaffected member is.
        ``dominant`` means that every affected member carries the variant and no unaffected member does.
        ``present`` covers the rest.
        Each event excludes the ones before it, so ``dominant`` means inherited dominant.
        In a trio with unaffected parents, ``dominant`` is therefore empty.
        Members of unknown status constrain nothing.
        Without an affected member, ``recessive`` and ``dominant`` are the impossible event ``false``.
        A pedigree without a child has no de novo event.
        One FDR control over ``present``, ``denovo_*``, ``recessive`` and ``dominant`` therefore fits every pedigree and keeps every variant.
        One over ``denovo_*`` and ``dominant`` fits the dominant hypothesis.
        The call passes ``--full-prior`` (``full_prior``), so that the mendelian prior and ``germline_mutation_rate`` weigh the events.
        Varlociraptor's default prior only tells absent from present, and it gives a child with unobserved parents a fixed share of ``denovo``.
        A member whose sequenced parent is outside the task is an error.
        This happens when ``library_selection`` drops the parent.
        Take both parents into the task, or give a ``template``.
        Libraries without a pedigree (``group_by: library``) are called alone.

    ``somatic``
        The scenario holds one tumor with its matched normal (``tumor_normal_mode``, as for ``mutect2``).
        It is the default scenario of dna-seq-varlociraptor.
        Its events are ``somatic_tumor_low``, ``somatic_tumor_medium`` and ``somatic_tumor_high`` (tumor allele frequency below 0.1, from 0.1 to 0.3, from 0.3), ``somatic_normal``, ``germline_het`` and ``germline_hom``.
        The tumor purity comes from the sample sheet column ``purity``, from 0 to 1.
        The scenario models it as the contamination of the tumor by the normal.
        A tumor without a purity counts as pure, with a warning.
        A tumor without a matched normal keeps the normal in the scenario as an unobserved sample named ``normal``, so the events stay the same.
        The VCF then has a sample ``normal`` without observations.
        The preset calls one tumor with its normal, so a donor with several tumors in the task (a primary and a relapse) is an error.
        Call each time point in a task of its own (``library_selection``), or give a ``template`` with ``group_by: related``.
        The template then receives the donor's tumors and the normal, and the entity is the primary tumor.
        This allows a joint scenario such as the tumor/normal/relapse one of the varlociraptor scenario catalog.

    In both presets, an FFPE library gets varlociraptor's model of the artifacts of an FFPE fragment (``fragment_end_bias``).
    The sample sheet column ``is_ffpe`` marks FFPE libraries.
    The samples of the scenario are the library names, and so are the samples of the VCF.
    The VCF carries ``GT`` from ``varlociraptor genotype``.
    The calls are all candidates with the posterior probability of each event (``INFO/PROB_<EVENT>``), not a filtered call set.
    Chain a ``variant_filtration`` task with ``tool: varlociraptor`` before any consumer, for FDR control on the events of interest.
    ``preprocess.max_depth`` caps the observations per library and position at 200 by downsampling.
    Raise it where the expected allele frequencies are low, such as for somatic variants of an impure tumor or for mosaics.
    A higher cap costs runtime in preprocessing and calling.
    A ``template`` is rendered with the variables ``samples``, ``library_name`` (the entity) and ``settings`` (``species`` and ``purity``).
    ``samples`` holds the library table rows of the entity's libraries, as dicts with ``library_name``, ``role``, ``sex``, ``father_name``, ``mother_name``, ``disease_state``, ``matched_normal_lib`` and the sheet's other columns.
    For the germline scenario, ``samples`` also holds the unsequenced relatives, with ``observed: false``.
    The sample keys of the template must be the library names.

By default (``group_by: related``), the germline callers call each pedigree jointly and produce one result per pedigree.
The result is named after the pedigree, its ``related_group``: the family id of the sample sheet (``pedigree_field``, else ``familyId``), else ``FAM_`` and the index donor.
With ``group_by: library``, they call each library on its own and produce one result per library.
``mutect2`` produces one result per tumor DNA library.
``library_selection`` can only narrow these.
``varlociraptor`` groups like the germline callers with ``scenario: germline`` and like ``mutect2`` with ``scenario: somatic``.
The GATK callers need ``depends_on.dbsnp`` in addition to ``depends_on.reference``.
With ``mutect2``, the optional ``depends_on.panel_of_normals`` names the panel VCF that Mutect2 uses to filter recurrent artefacts.
The panel comes from a ``panel_of_normals`` task with ``tool: mutect2``, or from an ``external_data`` task with ``{type: models, tags: [pon, mutect2]}`` and the key ``vcf``.
An example is ``1000g_pon.hg38.vcf.gz`` of the GATK resource bundle.
The optional ``depends_on.germline_resource`` names an ``external_data`` task with the population allele frequencies (``{type: variants, tags: [population]}``, key ``vcf``).
An example is the af-only gnomAD VCF of the GATK resource bundle.
With ``mutect2.contamination.enabled``, ``depends_on.common_variants`` is required: an ``external_data`` task with ``{type: variants, tags: [common]}`` and the key ``vcf``.
An example is ``small_exac_common_3.vcf`` of the GATK resource bundle.
Windows over the genome (GATK 4 ``PreprocessIntervals``) are processed in parallel.
The GATK callers and ``bcftools_call`` take ``window_length``, and these and ``mutect2`` take ``padding``, the bases by which neighbouring windows overlap.
The sections of these callers and of ``freebayes`` have ``ignore_contigs``, regular expressions for contigs that they skip.
By default, these are decoys, HLA, EBV and the hs37d5 extras.
A YAML anchor in ``config.yaml`` shares one list between tasks.
The ``FILTER`` IDs of the callers are reserved: ``LowQual`` of the GATK callers, and ``germline``, ``contamination`` and the other filters of ``FilterMutectCalls``.
No other task may be named like one of them, nor like ``PASS``.

Outputs
-------

Below ``tasks/<task>/output/``, for each ``{library_name}``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.
  For ``mutect2``, this file holds the filtered (PASS) calls.
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
    One job per task also summarizes the statistics of all samples in the result file ``qc_summary/out/qc_summary.tsv``.
    The summary has one row per sample and variant set, with genotype counts, het/hom and Ts/Tv ratios and depth.
    The job also renders tables and charts for the report below ``work/qc_summary/report/html/``: substitution types, indel length, depth and quality distributions.

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
      - name: reads
        step: external_data
        config:
          produces: {type: raw, tags: [fastq]}
          library_files: [reads.tsv]
      - name: mapping
        step: ngs_mapping
        config:
          depends_on: {reads: reads, reference: genome, index: bwa_index}
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

Varlociraptor on tumor/normal pairs of an exome cohort: candidates from freebayes in the kit's targets, calls against the somatic scenario, then FDR control on the somatic events.
The annotation and the other consumers read the result of the FDR control:

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
      - name: reads
        step: external_data
        config:
          produces: {type: raw, tags: [fastq]}
          library_files: [reads.tsv]
      - name: mapping
        step: ngs_mapping
        config:
          depends_on: {reads: reads, reference: genome, index: bwa_index}
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

For a germline cohort, ``candidates`` gets ``group_by: related``, ``calling`` gets ``scenario: germline`` and ``fdr_control`` the events ``present``, ``denovo_*`` (every ``denovo_<member>``), ``recessive`` and ``dominant``.
"""

import math
from pathlib import Path
from typing import Any

from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    MATCHED_NORMAL_RELATIONSHIP,
    MATCHED_NORMAL_REQUIREMENT,
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.common import vcf_outputs
from snappy_pipeline.workflows.abstract.common import (
    SnakemakeDict,
    SnakemakeDictItemsGenerator,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.common.scatter import RegionScatter, ScatterGather
from snappy_pipeline.workflows.variant_calling.model import ScenarioKind, Tool, TumorNormalMode
from snappy_pipeline.workflows.variant_calling.scenario import (
    NO_PARENT,
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


#: The germline callers call the normal DNA libraries
GERMLINE_REQUIREMENT = "role == 'normal' and extraction_type == 'dna'"


class VariantCallingStepPart(BaseStepPart):
    """Base class of the germline callers, which call the libraries of an entity jointly: the
    pedigree members of its related group with ``group_by: related``, else the library."""

    #: Class available actions
    actions = ("run",)

    result_keys = ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")

    def _get_log_file(self, action: str) -> SnakemakeDict:
        return {"log": f"work/{{library_name}}/log/{{library_name}}.{action}.log"}

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        rows = self._called_rows(wildcards.library_name)
        yield "ped", self._ped(wildcards.library_name, rows)
        yield (
            "bam",
            [
                self.parent.get_upstream_paths("alignments", library_name=name).bam
                for name in rows["library_name"]
            ],
        )

    def _get_output_files_run(self) -> dict[str, str]:
        return vcf_outputs("work/{library_name}/out/{library_name}")

    def members(self, library_name: str) -> list[str]:
        """Return the libraries called for the entity ``library_name``: the pedigree members of its related group with ``group_by: related``, else the library."""
        return self.parent.member_libraries(library_name)

    def _called_rows(self, library_name: str):
        """The library table rows of the libraries called for the entity, in table order"""
        df = self.parent.build_library_dataframe()
        return df[df["library_name"].isin(self.members(library_name))]

    def _ped(self, library_name: str, rows) -> str | list:
        """The PED file of the entity's family; none for libraries without one (cancer)"""
        if rows["family_id"].isna().all():
            return []
        if self.parent.effective_group_by != "library":
            group = library_name
        else:
            group = rows.iloc[0]["family_id"]
        return f"work/write_pedigree.{group}/out/{group}.ped"


class BcftoolsCallStepPart(VariantCallingStepPart):
    """Germline variant calling with bcftools"""

    #: Step name
    name = "bcftools_call"

    wrappers = {"run": "bcftools_call"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config.bcftools_call.num_threads
        mem_mb = int(3.75 * 1024 * num_threads)
        return ResourceUsage(threads=num_threads, runtime="2d", mem=f"{mem_mb}MB")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        """The options of the section; ``assembly`` names the ploidy of ``bcftools call``"""
        return self.config.bcftools_call.model_dump(mode="json")


class GatkCallerStepPartBase(VariantCallingStepPart):
    """Base class for GATK v3/v4 variant callers"""

    @property
    def default_resource_usage(self) -> ResourceUsage:
        num_threads = self.config[self.name].num_threads
        mem_per_thread = 5.5
        mem_total = int(mem_per_thread * num_threads + 0.5)
        return ResourceUsage(threads=num_threads, runtime="2d", mem=f"{mem_total}GB")

    @dictify
    def _get_input_files_run(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        """The BAMs and the pedigree, the reference and dbSNP"""
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        yield from self.parent.vcf_inputs("dbsnp").items()


class GatkJointCallerStepPart(GatkCallerStepPartBase):
    """Base class of the GATK callers that call an entity in one job; the params are the options
    of the tool's section"""

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return self.tool_section().model_dump()


class Gatk3HaplotypeCallerStepPart(GatkJointCallerStepPart):
    """Germline variant calling with GATK v3 HaplotypeCaller"""

    #: Step name
    name = "gatk3_hc"

    wrappers = {"run": "gatk3_hc"}


class Gatk3UnifiedGenotyperStepPart(GatkJointCallerStepPart):
    """Germline variant calling with GATK v3 UnifiedGenotyper"""

    #: Step name
    name = "gatk3_ug"

    wrappers = {"run": "gatk3_ug"}


class Gatk4HaplotypeCallerJointStepPart(GatkJointCallerStepPart):
    """Germline variant calling with GATK 4 HaplotypeCaller doing joint calling per pedigree"""

    name = "gatk4_hc_joint"
    wrappers = {"run": "gatk4_hc/joint"}


class Gatk4HaplotypeCallerGvcfStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller and gVCF workflow"""

    name = "gatk4_hc_gvcf"
    actions = ("discover", "combine_gvcfs", "genotype")

    wrappers = {
        "genotype": "gatk4_hc/genotype",
        "combine_gvcfs": "gatk4_hc/combine_gvcfs",
        "discover": "gatk4_hc/discover",
    }

    extra_arguments_action = "discover"

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

    intermediate_outputs = {"discover": ("gvcf", "gvcf_tbi")}

    #: The gVCF of a library, and the combined gVCF of an entity
    discovered_tpl = "work/{library_name}/out/{library_name}.discover.g.vcf.gz"

    combined_tpl = "work/{library_name}/out/{library_name}.combine_gvcfs.g.vcf.gz"

    @dictify
    def _get_input_files_discover(self, wildcards):
        yield from self.parent.reference_inputs().items()
        yield from self.parent.vcf_inputs("dbsnp").items()
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam

    def _get_output_files_discover(self) -> dict[str, str]:
        return {"gvcf": self.discovered_tpl, "gvcf_tbi": self.discovered_tpl + ".tbi"}

    def _get_params_discover(self, wildcards) -> dict[str, Any]:
        return self._params_gvcf("discover")

    @dictify
    def _get_input_files_combine_gvcfs(self, wildcards: Wildcards) -> SnakemakeDictItemsGenerator:
        yield from self.parent.reference_inputs().items()
        rows = self._called_rows(wildcards.library_name)
        yield "ped", self._ped(wildcards.library_name, rows)
        gvcfs = [self.discovered_tpl.format(library_name=name) for name in rows["library_name"]]
        yield "gvcf", gvcfs
        yield "gvcf_tbi", [f"{gvcf}.tbi" for gvcf in gvcfs]

    def _get_output_files_combine_gvcfs(self) -> dict[str, str]:
        return {"gvcf": self.combined_tpl, "gvcf_tbi": self.combined_tpl + ".tbi"}

    def _get_params_combine_gvcfs(self, wildcards) -> dict[str, Any]:
        return self._params_gvcf("combine_gvcfs")

    @dictify
    def _get_input_files_genotype(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield from self.parent.reference_inputs().items()
        gvcf = self.combined_tpl.format(**wildcards)
        yield "gvcf", gvcf
        yield "gvcf_tbi", gvcf + ".tbi"

    def _get_output_files_genotype(self) -> dict[str, str]:
        prefix = "work/{library_name}/out/{library_name}"
        gvcf = {f"g{key}": path for key, path in vcf_outputs(f"{prefix}.g").items()}
        return gvcf | vcf_outputs(prefix)

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


class Mutect2StepPart(RegionScatter, BaseStepPart):
    """Somatic calling with Mutect2: each tumor library with its matched normal, scattered over
    the regions of the genome, which the task's tumors share"""

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

    @property
    def wrappers(self) -> dict[str, str]:
        """The contamination estimate and its pileups run only when ``contamination`` is enabled"""
        wrappers = RegionScatter.scatter_wrapper | {
            "filter": "mutect2/filter",
            "gather": "mutect2/gather",
            "run": "mutect2/run",
        }
        if self.config.mutect2.contamination.enabled:
            wrappers |= {
                "contamination": "mutect2/contamination",
                "pileup_tumor": "mutect2/pileup",
                "pileup_normal": "mutect2/pileup",
            }
        return wrappers

    #: The sub-sections ``filtration`` and ``contamination`` hold the arguments of the other
    #: GATK actions
    extra_arguments_action = "run"

    result_keys = {"filter": ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")}
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
            "run": ResourceUsage(threads=self.config.mutect2.num_cores, runtime="5d", mem="8000MB"),
            "gather": ResourceUsage(threads=1, runtime="4h", mem="32768MB"),
            "filter": ResourceUsage(threads=2, runtime="4h", mem="15872MB"),
            "contamination": ResourceUsage(threads=2, runtime="4h", mem="7680MB"),
            "pileup_normal": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
            "pileup_tumor": ResourceUsage(threads=2, runtime="4h", mem="8000MB"),
        }

    @property
    def scatter_padding(self) -> int:
        return self.config.mutect2.padding

    #: The raw calls of a tumor in one part of the genome
    chunk_tpl = "work/{library_name}/out/{library_name}/calls/{scatteritem}.raw"

    #: The outputs of a part and their extensions
    _run_outputs = {
        "vcf": ".vcf.gz",
        "vcf_tbi": ".vcf.gz.tbi",
        "stats": ".vcf.stats",
        "f1r2": ".f1r2.tar.gz",
    }

    def _get_log_file(self, action):
        if action == "run":
            postfix = ".{scatteritem}"
        elif action == "gather":
            postfix = ""
        else:
            postfix = "." + action
        return {"log": f"work/{{library_name}}/log/{{library_name}}{postfix}.log"}

    def scatter_ignore_contigs(self) -> list[str]:
        return sorted(set(self.config.mutect2.ignore_contigs))

    def _get_input_files_run(self, wildcards):
        tumor = self.parent.get_upstream_paths("alignments", library_name=wildcards.library_name)
        input_files = {
            "tumor_bam": tumor.bam,
            "tumor_bai": tumor.bai,
            "region": self.region(wildcards),
        }
        # ``paired`` mode takes only the tumors with a matched normal (task_library_requirement)
        normal_library = self.get_normal_lib_name(wildcards)
        if normal_library and self.config.mutect2.tumor_normal_mode != TumorNormalMode.TUMOR_ONLY:
            normal = self.parent.get_upstream_paths("alignments", library_name=normal_library)
            input_files.update({"normal_bam": normal.bam, "normal_bai": normal.bai})

        input_files.update(self.parent.reference_inputs())

        if self.config.depends_on.germline_resource:
            input_files.update(self.parent.vcf_inputs("germline_resource"))
        if self.config.depends_on.panel_of_normals:
            input_files.update(self.parent.vcf_inputs("panel_of_normals"))

        return input_files

    def _get_output_files_run(self):
        return {key: self.chunk_tpl + ext for key, ext in self._run_outputs.items()}

    def _get_params_run(self, wildcards):
        return {
            "normal_lib_name": self.get_normal_lib_name(wildcards),
            "java_options": self.config.mutect2.java_options,
            "extra_arguments": self.config.mutect2.extra_arguments,
        }

    def _get_input_files_gather(self, wildcards):
        return {
            key: self.gathered(self.chunk_tpl + ext, wildcards)
            for key, ext in self._run_outputs.items()
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

    def _get_output_files_filter(self):
        return vcf_outputs("work/{library_name}/out/{library_name}") | self._with_exts(
            {"full_vcf": ".full.vcf.gz", "full_vcf_tbi": ".full.vcf.gz.tbi"}
        )

    def _get_params_filter(self, wildcards):
        return self.config.mutect2.filtration.model_dump(by_alias=True)

    def _get_input_files_contamination(self, wildcards: Wildcards):
        base_path = "work/{library_name}/out/{library_name}".format(**wildcards)
        return {"normal": base_path + ".normal.pileup", "tumor": base_path + ".tumor.pileup"}

    def _get_output_files_contamination(self):
        return self._with_exts({"table": ".contamination.tbl", "segments": ".segments.tbl"})

    def _get_params_contamination(self, wildcards):
        return {
            "java_options": self.config.mutect2.contamination.java_options,
            "extra_arguments": self.config.mutect2.contamination.extra_arguments,
        }

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

    def _get_output_files_pileup_normal(self):
        return self._with_exts({"pileup": ".normal.pileup"})

    def _get_params_pileup_normal(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "normal_lib_name": self.get_normal_lib_name(wildcards)
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

    def _get_output_files_pileup_tumor(self):
        return self._with_exts({"pileup": ".tumor.pileup"})

    def _get_params_pileup_tumor(self, wildcards):
        return self.config.mutect2.contamination.pileup.model_dump(by_alias=True) | {
            "tumor_lib_name": self.get_tumor_lib_name(wildcards)
        }

    def get_normal_lib_name(self, wildcards):
        """Return the matched normal library of the tumor library, or None"""
        return self.parent.matched_normal(wildcards.library_name) or None

    def members(self, library_name: str) -> list[str]:
        """Return the libraries in the VCF of the tumor ``library_name``: the tumor and its matched normal, if it has one."""
        normal = self.parent.matched_normal(library_name)
        return [library_name, normal] if normal else [library_name]

    def get_tumor_lib_name(self, wildcards):
        return wildcards.library_name

    @staticmethod
    def _with_exts(exts: dict[str, str]) -> dict[str, str]:
        return {key: "work/{library_name}/out/{library_name}" + ext for key, ext in exts.items()}


class ScatteredCallerStepPart(ScatterGather, BaseStepPart):
    """Base class of the callers that call an entity in scattered parts and concatenate them
    with the ``bcftools/concat`` wrapper."""

    #: The calls of one part of an entity, a BCF with its ``.csi``
    chunk_tpl = "work/{library_name}/out/{library_name}/calls/{scatteritem}"

    def _get_input_files_gather(self, wildcards):
        return {
            "vcf": self.gathered(self.chunk_tpl + ".bcf", wildcards),
            "vcf_csi": self.gathered(self.chunk_tpl + ".bcf.csi", wildcards),
        }

    def _get_output_files_gather(self):
        return vcf_outputs("work/{library_name}/out/{library_name}")

    def _get_log_file_gather(self):
        return {"log": "work/{library_name}/log/{library_name}.gather.log"}

    def members(self, library_name: str) -> list[str]:
        """Return the libraries called together for the entity ``library_name``: with ``group_by: related`` the pedigree members of its related group (one DNA library per donor, the index first), else the library and its matched normal, if it has one."""
        if self.parent.effective_group_by != "library":
            return self.parent.member_libraries(library_name)
        normal = self.parent.matched_normal(library_name)
        return [library_name, normal] if normal else [library_name]


class FreebayesStepPart(RegionScatter, ScatteredCallerStepPart):
    """Candidate calling with freebayes, jointly over the libraries of an entity, per genome
    window; the windows are shared by the task's entities and disjoint (no padding), so no call
    is made twice"""

    name = "freebayes"
    actions = ("scatter", "run", "gather")
    wrappers = RegionScatter.scatter_wrapper | {"gather": "bcftools/concat", "run": "freebayes"}
    extra_arguments_action = "run"
    result_keys = {"gather": ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")}
    temp_outputs = {"scatter": ("regions",), "run": ("bcf", "bcf_csi")}

    resource_usage = {
        "run": ResourceUsage(threads=1, runtime="1d", mem="8GB"),
        "gather": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
    }

    def scatter_ignore_contigs(self) -> list[str]:
        return sorted(set(self.config.freebayes.ignore_contigs))

    def _get_input_files_scatter(self, wildcards):
        """The genome index and, with ``depends_on.kits``, the targets of the kits of the
        task's libraries: the windows then hold equal target bases instead of equal lengths"""
        inputs = super()._get_input_files_scatter(wildcards)
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
        yield "region", self.region(wildcards)
        if self.config.depends_on.kits:
            kit = self.parent.common_kit("kits", members)
            yield "targets", self.parent.get_upstream_paths("kits", library_kit=kit).targets

    def _get_output_files_run(self):
        return {"bcf": self.chunk_tpl + ".bcf", "bcf_csi": self.chunk_tpl + ".bcf.csi"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.run.{scatteritem}.log"}

    def _get_params_run(self, wildcards):
        return {
            "min_alternate_fraction": self.config.freebayes.min_alternate_fraction,
            "min_alternate_count": self.config.freebayes.min_alternate_count,
            "atomize": self.config.freebayes.atomize,
        }


class VarlociraptorStepPart(ScatteredCallerStepPart):
    """Calling with varlociraptor on the candidates of ``depends_on.candidates``: the candidates
    of an entity split into chunks of equal depth, alignment properties per library,
    observations per library and chunk, and calls per chunk against the scenario of the entity"""

    name = "varlociraptor"
    actions = ("scatter", "alignment_properties", "preprocess", "call", "gather")

    wrappers = {
        "gather": "bcftools/concat",
        "call": "varlociraptor/call",
        "preprocess": "varlociraptor/preprocess",
        "alignment_properties": "varlociraptor/alignment_properties",
        "scatter": "varlociraptor/split_candidates",
    }

    #: The ``preprocess`` section passes its own extra arguments
    extra_arguments_action = "call"

    result_keys = {"gather": ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")}

    temp_outputs = {
        "scatter": ("candidates",),
        "preprocess": ("bcf",),
        "call": ("bcf", "bcf_csi"),
    }

    intermediate_outputs = {"alignment_properties": ("json",)}

    resource_usage = {
        "scatter": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
        "alignment_properties": ResourceUsage(threads=1, runtime="2h", mem="4GB"),
        # The joint model grows with the samples of the pedigree.
        # A chunk of an eight-member family with complex candidates took 14 GB.
        "preprocess": ResourceUsage(threads=1, runtime="1d", mem="8GB"),
        "call": ResourceUsage(threads=1, runtime="1d", mem="16GB"),
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
        tumor = self.parent.primary_library(library_name)
        if cfg.scenario == ScenarioKind.germline:
            if self.parent.effective_group_by != "library":
                return super().members(library_name)
            return [tumor]
        members = [tumor]
        if self.parent.effective_group_by != "library":
            members += self.other_tumors(library_name)
        if cfg.tumor_normal_mode != TumorNormalMode.TUMOR_ONLY:
            if normal := self.parent.matched_normal(tumor):
                members.append(normal)
        return members

    def _get_input_files_scatter(self, wildcards):
        candidates = self.parent.get_upstream_paths(
            "candidates", library_name=wildcards.library_name
        )
        return {"candidates": candidates.vcf}

    def _get_output_files_scatter(self):
        template = self.candidates_tpl.replace("{library_name}", "{{library_name}}")
        return {"candidates": self._scattergather("scatter")(template)}

    def _get_log_file_scatter(self):
        return {"log": "work/{library_name}/log/{library_name}.scatter.log"}

    def _get_input_files_alignment_properties(self, wildcards):
        alignments = self.parent.get_upstream_paths("alignments", library_name=wildcards.member)
        return {**self.parent.reference_inputs(), "bam": alignments.bam, "bai": alignments.bai}

    def _get_output_files_alignment_properties(self):
        return {"json": self.properties_tpl}

    def _get_log_file_alignment_properties(self):
        return {"log": "work/{member}/log/{member}.alignment_properties.log"}

    def _get_input_files_preprocess(self, wildcards):
        alignments = self.parent.get_upstream_paths("alignments", library_name=wildcards.member)
        return {
            **self.parent.reference_inputs(),
            "candidates": self.candidates_tpl.format(**wildcards),
            "alignment_properties": self.properties_tpl.format(**wildcards),
            "bam": alignments.bam,
            "bai": alignments.bai,
        }

    def _get_output_files_preprocess(self):
        return {"bcf": self.observations_tpl}

    def _get_log_file_preprocess(self):
        return {
            "log": "work/{library_name}/log/{library_name}.preprocess.{member}.{scatteritem}.log"
        }

    def _get_params_preprocess(self, wildcards):
        preprocess = self.config.varlociraptor.preprocess
        return {"max_depth": preprocess.max_depth, "extra_arguments": preprocess.extra_arguments}

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

    def _get_output_files_call(self):
        return {"bcf": self.chunk_tpl + ".bcf", "bcf_csi": self.chunk_tpl + ".bcf.csi"}

    def _get_log_file_call(self):
        return {"log": "work/{library_name}/log/{library_name}.call.{scatteritem}.log"}

    def _get_params_call(self, wildcards):
        """The sample names of the observations in their order, the scenario,
        and whether varlociraptor uses the full prior"""
        return {
            "obs_names": self.members(wildcards.library_name),
            "scenario": self.scenario(wildcards.library_name),
            "full_prior": self.config.varlociraptor.full_prior,
        }

    def other_tumors(self, library_name: str) -> list[str]:
        """The other tumor libraries of the donor of the entity ``library_name`` that the task takes"""
        libraries = self.parent.libraries
        tumor = self.parent.primary_library(library_name)
        same_donor = libraries["donor_name"] == self.parent.library_row(tumor)["donor_name"]
        other = libraries["library_name"] != tumor
        return libraries[same_donor & other]["library_name"].tolist()

    def sample(self, library_name: str) -> dict[str, Any]:
        """The library table row of a library with builtin values, as a scenario template sees it"""
        return {
            key: _builtin(value) for key, value in self.parent.library_row(library_name).items()
        }

    def pedigree_samples(self, library_name: str) -> list[dict[str, Any]]:
        """The samples of the germline scenario of a pedigree.
        The observed members come first.
        Unsequenced donors that connect members join as samples without observations,
        named after the donor: a parent of two sequenced children, the grandparents of cousins.
        A sample with one parent in the scenario and an unknown other parent gets an unobserved parent.
        It is named ``<donor>.father`` or ``<donor>.mother``, since varlociraptor inherits from both parents or none.
        Parents are the sample names, or ``"0"``.
        A sequenced parent outside the task keeps its library name, so that the preset rejects it."""
        members = self.members(library_name)
        samples = [self.sample(name) for name in members]
        if self.parent.effective_group_by == "library":
            return samples
        donors = self.parent.pedigree_donors(library_name, include_background=True)
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
                sample.setdefault("donor_name", donor)  # names the member's de novo event
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
        return samples + self.unknown_parents(samples, by_library)

    @staticmethod
    def unknown_parents(
        samples: list[dict[str, Any]], by_library: dict[str, str]
    ) -> list[dict[str, Any]]:
        """Return an unobserved parent for every sample with one parent among ``samples`` and an unknown other parent.
        It is named ``<donor>.father`` or ``<donor>.mother`` and set as that parent.
        A parent outside the samples, a sequenced one that the task leaves out, stays for the preset to reject."""
        names = {sample["library_name"] for sample in samples}
        parents = []
        for sample in samples:
            father, mother = sample["father_name"], sample["mother_name"]
            if (father in names) == (mother in names) or not {father, mother} & set(NO_PARENT):
                continue
            key, sex = ("mother_name", "female") if father in names else ("father_name", "male")
            donor = by_library.get(sample["library_name"], sample["library_name"])
            name = f"{donor}.{key.removesuffix('_name')}"
            if name in names:
                raise ValueError(f"{name} is a sample and the unknown {key[:-5]} of {donor}")
            sample[key] = name
            parents.append(
                {
                    "library_name": name,
                    "sex": sex,
                    "father_name": "0",
                    "mother_name": "0",
                    "observed": False,
                }
            )
        return parents

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
        purity = self.purity(members[0]) if somatic else None
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
                    f"of {members[0]} has the tumors {', '.join([members[0], *others])} in "
                    "the task; call each time point in a task of its own (library_selection), or "
                    "give a scenario template with group_by: related for a joint scenario of the "
                    "donor's tumors"
                )
            normal = members[1] if len(members) > 1 else None
            ffpe = frozenset(name for name in members if is_ffpe(self.sample(name).get("is_ffpe")))
            scenario = somatic_scenario(members[0], normal, purity, ffpe)
        else:
            scenario = germline_scenario(
                self.pedigree_samples(library_name), cfg.species.model_dump()
            )
        check_scenario(scenario, members)
        return scenario


class ReportStepPart(BaseStepPart):
    """Base class for the reports on the called variants, one per entity (and donor)"""

    actions = ("run",)
    entity_wildcard = "library_name"

    #: Whether we generate per-donor files.
    report_per_donor: bool

    def _get_log_file(self, action: str) -> SnakemakeDict:
        token = "{library_name}"
        donor = "{donor_library_name}." if self.report_per_donor else ""
        return {"log": f"work/{token}/log/{token}.{donor}{self.name}_{action}.log"}

    def _result_wildcards(self, entity: str) -> dict[str, list[str]]:
        result = super()._result_wildcards(entity)
        if self.report_per_donor and self.parent.effective_group_by != "library":
            result["donor_library_name"] = self.parent.called_libraries(entity)
        return result


class BafFileGenerationStepPart(ReportStepPart):
    """Class for computing B allele fraction files.

    One file is generated per sample in the output VCF files.
    """

    #: Step name
    name = "baf_file_generation"

    wrappers = {"run": "baf_file_generation"}
    result_keys = ("bw", "bw_md5")
    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="1024MB")

    #: Whether to report results per donor
    report_per_donor = True

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{library_name}/out/{library_name}.vcf.gz".format(**wildcards),
        )
        yield "reference_index", self.parent.get_upstream_paths("reference").fai

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        base_path = "{library_name}/report/baf/{library_name}.{donor_library_name,[^\\.]+}.baf"
        ext_names = {"bw": ".bw", "bw_md5": ".bw.md5"}
        work_files = {}
        for key, ext in ext_names.items():
            work_files[key] = f"work/{base_path}{ext}"
        yield from work_files.items()

    def _get_params_run(self, wildcards):
        return {"min_dp": self.config.baf_file_generation.min_dp}


class BcftoolsStatsStepPart(ReportStepPart):
    """VCF statistics computation with ``bcftools stats``.

    Statistics are computed overall and per-sample
    """

    name = "bcftools_stats"
    wrappers = {"run": "bcftools_stats/run"}
    result_keys = ("txt", "txt_md5")
    report_outputs = {"run": {"txt": ReportOutput("bcftools stats", caption="bcftools_stats.rst")}}
    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="1024MB")
    report_per_donor = True

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{library_name}/out/{library_name}.vcf.gz".format(**wildcards),
        )

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = "work/{library_name}/report/bcftools_stats/{library_name}.{donor_library_name}"
        work_files = {key: f"{base_path}{ext}" for key, ext in ext_names.items()}
        yield from work_files.items()

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {"donor_library_name": wildcards.donor_library_name}


class BcftoolsRohStepPart(ReportStepPart):
    """ROH calling with `bcftools roh`."""

    name = "bcftools_roh"
    wrappers = {"run": "bcftools_roh"}
    result_keys = ("txt", "txt_md5")
    default_resource_usage = ResourceUsage(threads=1, runtime="1h", mem="4GB")
    report_per_donor = False

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        yield (
            "vcf",
            "work/{library_name}/out/{library_name}.vcf.gz".format(**wildcards),
        )
        yield "path_af_file", self.config.get(self.name).get("path_af_file")
        if self.parent.config.depends_on.kits:
            kit = self.parent.common_kit(
                "kits", self.parent.called_libraries(wildcards.library_name)
            )
            yield "path_targets", self.parent.get_upstream_paths("kits", library_kit=kit).targets

    @dictify
    def _get_output_files_run(self) -> SnakemakeDictItemsGenerator:
        ext_names = {"txt": ".txt", "txt_md5": ".txt.md5"}
        base_path = "work/{library_name}/report/roh/{library_name}"
        for key, ext in ext_names.items():
            yield key, f"{base_path}{ext}"

    def _get_params_run(self, wildcards):
        return {
            name: self.config.bcftools_roh.get(name)
            for name in [
                "ignore_homref",
                "skip_indels",
                "rec_rate",
            ]
        }


class BcftoolsStatsSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the ``bcftools stats`` of all samples"""

    report_outputs = {
        "run": {
            "html": ReportOutput(
                "bcftools stats", caption="bcftools_stats_summary.rst", htmlindex="index.html"
            )
        }
    }

    source = ("bcftools_stats", "run")
    source_keys = {"txt": "txt"}
    format = "bcftools_stats"

    #: One file per sample of the entity's VCF
    label_wildcard = "donor_library_name"


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

    tool_parts = {
        Tool.bcftools_call: BcftoolsCallStepPart,
        Tool.gatk3_hc: Gatk3HaplotypeCallerStepPart,
        Tool.gatk3_ug: Gatk3UnifiedGenotyperStepPart,
        Tool.gatk4_hc_joint: Gatk4HaplotypeCallerJointStepPart,
        Tool.gatk4_hc_gvcf: Gatk4HaplotypeCallerGvcfStepPart,
        Tool.mutect2: Mutect2StepPart,
        Tool.freebayes: FreebayesStepPart,
        Tool.varlociraptor: VarlociraptorStepPart,
    }

    helper_parts = (WritePedigreeStepPart,)

    #: Default relationship for somatic (mutect2) callers: resolve matched normal.
    default_relationships = MATCHED_NORMAL_RELATIONSHIP

    #: The reports on the called variants, in the order of their result files
    report_parts = (BafFileGenerationStepPart, BcftoolsStatsStepPart, BcftoolsRohStepPart)

    @classmethod
    def task_library_requirement(cls, config):
        """mutect2 and the somatic scenario of varlociraptor call tumor DNA libraries, which
        need a matched normal or none; ``tumor_normal_mode: paired`` leaves out the tumors
        without one. The germline callers call the normal DNA libraries; freebayes finds the
        candidates of every DNA library."""
        if config.tool == Tool.mutect2:
            mode = config.mutect2.tumor_normal_mode
        elif config.tool == Tool.varlociraptor:
            if config.varlociraptor.scenario != ScenarioKind.somatic:
                return GERMLINE_REQUIREMENT
            mode = config.varlociraptor.tumor_normal_mode
        elif config.tool == Tool.freebayes:
            return None  # candidates of every DNA library, tumors included
        else:
            return GERMLINE_REQUIREMENT
        if mode == TumorNormalMode.PAIRED:
            return MATCHED_NORMAL_REQUIREMENT
        return "role == 'tumor' and extraction_type == 'dna'"

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

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        enabled = [part for part in self.report_parts if self._enabled(part.name)]
        if self._enabled("bcftools_stats"):
            enabled.append(BcftoolsStatsSummaryStepPart)
        self.register_sub_step_classes(tuple(enabled))
        if self.config.depends_on.kits and "bcftools_roh" in self.sub_steps:
            for library_name in self.output_entities:
                self.common_kit("kits", self.called_libraries(library_name))
        if self.config.tool == Tool.freebayes and self.config.depends_on.kits:
            # The candidates of an entity are restricted to the targets of its libraries' kit
            for library_name in self.output_entities:
                self.common_kit("kits", self.sub_steps["freebayes"].members(library_name))
        if self.config.tool == Tool.varlociraptor:
            self._check_varlociraptor_setup()

    def library_warnings(self) -> list[str]:
        """The tumors without a matched normal that ``paired`` mode leaves out, and the tumors
        that the somatic scenario of varlociraptor takes as pure for lack of a purity."""
        warnings = super().library_warnings()
        if (
            self.config.tool == Tool.varlociraptor
            and self.config.varlociraptor.scenario == ScenarioKind.somatic
        ):
            part = self.sub_steps["varlociraptor"]
            unknown = [
                name
                for name in self.output_entities
                if part.purity(self.primary_library(name)) is None
            ]
            if unknown:
                warnings.append(
                    f"Task {self.task_name!r} takes the tumors without a purity as pure, with no "
                    f"contamination by the normal: {', '.join(unknown)}; the sample sheet column "
                    "purity (a custom field, from 0 to 1) sets it, from the lab, the collaborator "
                    "or a purity estimate (somatic_purity_ploidy_estimate)"
                )
        return warnings

    def called_libraries(self, entity: str) -> list[str]:
        """Return the libraries in the VCF of an entity, as the task's caller takes them (its ``members``)."""
        return self.sub_steps[self.config.tool].members(entity)

    def _check_varlociraptor_setup(self) -> None:
        """The candidates task writes its calls per entity of this task, and the library names
        are sample names that a scenario accepts."""
        candidates = self.resolve_dependency("candidates").task_name
        upstream_group_by = self.upstream_group_by(candidates)
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

    def _enabled(self, report: str) -> bool:
        """Whether the config enables the report"""
        section = self.config.get(report)
        return bool(section and section.enabled)

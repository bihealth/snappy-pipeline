# -*- coding: utf-8 -*-
"""Small variant calling

Calls small variants (SNVs and indels) on the alignments of the ``ngs_mapping`` step.
The ``tool`` field selects the caller.

``gatk4_hc_gvcf``
    GATK 4 HaplotypeCaller in the gVCF workflow: variant discovery per library, ``CombineGVCFs`` per cohort, then ``GenotypeGVCFs``.
    The output VCF carries physical phasing information.
    The task also publishes the gVCF files.

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

The germline callers produce one result per cohort with ``group_by: cohort``, named after the index library, and one per library with ``group_by: library``.
``mutect2`` produces one result per tumor library.
The GATK callers need ``depends_on.dbsnp`` in addition to ``depends_on.reference``.
With ``mutect2``, ``depends_on.panel_of_normals`` optionally names a ``panel_of_normals`` task with ``tool: mutect2``, whose panel VCF Mutect2 uses to filter recurrent artefacts.
``depends_on.germline_resource`` optionally names an ``external_data`` task with the population allele frequencies (``{type: variants, tags: [population]}``, key ``vcf``, for example the af-only gnomAD VCF of the GATK resource bundle).
With ``mutect2.contamination.enabled``, ``depends_on.common_variants`` is required: an ``external_data`` task with ``{type: variants, tags: [common]}`` and the key ``vcf``, for example ``small_exac_common_3.vcf`` of the GATK resource bundle.
Windows over the genome (GATK 4 ``PreprocessIntervals``) are processed in parallel.

Outputs
-------

Below ``tasks/<task>/output/``, for each ``{library_name}``:

- ``{library_name}/out/{library_name}.vcf.gz`` with ``.tbi``, ``.md5`` and ``.tbi.md5`` files.
  For ``mutect2`` this is the filtered (PASS) calls.
- ``{library_name}/out/{library_name}.g.vcf.gz`` with the same companion files, for ``gatk4_hc_gvcf`` only.

Downstream tasks read ``vcf`` and ``vcf_tbi`` through ``depends_on``.
The germline callers produce ``variants`` with tags ``germline``, ``snv``, ``indel``.
``mutect2`` produces ``somatic``, ``snv``, ``indel`` for all calls and additionally ``filtered`` for the PASS calls.

Statistics
----------

The enabled statistics are published below ``{library_name}/report/``, with ``.md5`` files.

``bcftools_stats``
    ``bcftools stats`` per donor of each cohort: ``bcftools_stats/{library_name}.{donor_library_name}.txt``.

``bcftools_roh``
    Runs of homozygosity with ``bcftools roh``: ``roh/{library_name}.txt``.
    With ``depends_on.kits``, the regions are the ``targets`` of the kit of the targeted libraries in the VCF, which must share one kit.
    Without it, ROH are called on the whole genome.

``baf_file_generation``
    B-allele fraction per donor as a BigWig file, to confirm structural variants or runs of homozygosity visually.
"""

from typing import Any

from biomedsheets.shortcuts import GermlineCaseSheet
from snakemake.iocontainers import Wildcards

from snappy_pipeline.models import RelationshipDefinition
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ReportOutput,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import (
    SnakemakeDict,
    SnakemakeDictItemsGenerator,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.ngs_mapping.model import ExpectedAlignments
from snappy_pipeline.workflows.variant_calling.model import Tool, TumorNormalMode

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

#: Default configuration for the variant_calling step


class VariantCallingStepPart(BaseStepPart):
    """Base class for germline variant calling step parts

    Variant calling is performed on a per-pedigree level.  The (single) index individual is used
    for naming the output file.
    """

    #: Class available actions
    actions = ("run",)

    published_keys = ("vcf", "vcf_md5", "vcf_tbi", "vcf_tbi_md5")

    def _get_log_file(self, action: str) -> SnakemakeDict:
        token = "{library_name}"
        return {"log": f"work/{token}/log/{token}.{self.name}_{action}.log"}

    @dictify
    def _get_input_files_run(self, wildcards) -> SnakemakeDictItemsGenerator:
        df = self.parent.build_library_dataframe()
        if df.empty:
            return

        group_by = self.parent.effective_group_by or "cohort"

        # Determine the target libraries based on the wildcard
        if group_by == "cohort":
            # wildcards.library_name is the index library name; find its cohort
            df_target = df[
                df["library_name"].isin(self.parent.get_cohort_libraries(wildcards.library_name))
            ]
        elif group_by == "library":
            # Per-library logic
            df_target = df[df["library_name"] == wildcards.library_name]
        else:
            raise ValueError(f"Unknown group_by: {group_by}")

        selection = getattr(self.config, "library_selection", None)
        if selection:
            from snappy_pipeline.workflows.abstract import apply_library_selection

            # Note: "kind" logic isn't strictly necessary if selection encompasses it, but we can assume default 'germline' if needed
            df_target = apply_library_selection(df_target, selection, "germline")

        if df_target["kind"].eq("cancer").all():
            yield "ped", []
        else:
            if group_by == "cohort":
                yield (
                    "ped",
                    f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
                )
            else:
                cohort_name = df_target.iloc[0]["cohort_name"]
                yield "ped", f"work/write_pedigree.{cohort_name}/out/{cohort_name}.ped"

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
            "ignore_chroms": self.config.ignore_chroms,
            "window_length": self.config.bcftools_call.window_length,
            "num_threads": self.config.bcftools_call.num_threads,
            "max_depth": self.config.bcftools_call.max_depth,
            "max_indel_depth": self.config.bcftools_call.max_indel_depth,
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
        yield "dbsnp", self.parent.get_upstream_paths("dbsnp").vcf

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
            "ignore_chroms": self.config.ignore_chroms,
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
            "ignore_chroms": self.config.ignore_chroms,
        }


class Gatk4HaplotypeCallerJointStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller doing joint calling per pedigree"""

    name = "gatk4_hc_joint"

    def _get_params_run(self, wildcards) -> dict[str, Any]:
        return {
            "window_length": self.config.gatk4_hc_joint.window_length,
            "num_threads": self.config.gatk4_hc_joint.num_threads,
            "allow_seq_dict_incompatibility": self.config.gatk4_hc_joint.allow_seq_dict_incompatibility,
            "ignore_chroms": self.config.ignore_chroms,
        }


class Gatk4HaplotypeCallerGvcfStepPart(GatkCallerStepPartBase):
    """Germline variant calling with GATK 4 HaplotypeCaller and gVCF workflow"""

    name = "gatk4_hc_gvcf"

    actions = ("discover", "combine_gvcfs", "genotype")

    intermediate_outputs = {"discover": ("gvcf", "gvcf_tbi")}

    #: The genotype rule links the combined gVCF as ``{library}.g.vcf.gz``, for re-genotyping
    published_keys = {
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
        yield "dbsnp", self.parent.get_upstream_paths("dbsnp").vcf
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

        group_by = self.parent.effective_group_by or "cohort"

        # Determine the target libraries based on the wildcard
        if group_by == "cohort":
            # wildcards.library_name is the index library name; find its cohort
            df_target = df[
                df["library_name"].isin(self.parent.get_cohort_libraries(wildcards.library_name))
            ]
        elif group_by == "library":
            # Per-library logic
            df_target = df[df["library_name"] == wildcards.library_name]
        else:
            raise ValueError(f"Unknown group_by: {group_by}")

        selection = getattr(self.config, "library_selection", None)
        if selection:
            from snappy_pipeline.workflows.abstract import apply_library_selection

            df_target = apply_library_selection(df_target, selection, "germline")

        if df_target["kind"].eq("cancer").all():
            yield "ped", []
        else:
            if group_by == "cohort":
                yield (
                    "ped",
                    f"work/write_pedigree.{wildcards.library_name}/out/{wildcards.library_name}.ped",
                )
            else:
                cohort_name = df_target.iloc[0]["cohort_name"]
                yield "ped", f"work/write_pedigree.{cohort_name}/out/{cohort_name}.ped"

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
        return {
            "step_key": "variant_calling",
            "caller_key": "gatk4_hc_gvcf",
            "window_length": self.config.gatk4_hc_gvcf.window_length,
            "num_threads": self.config.gatk4_hc_gvcf.num_threads,
            "allow_seq_dict_incompatibility": self.config.gatk4_hc_gvcf.allow_seq_dict_incompatibility,
            "ignore_chroms": self.config.ignore_chroms,
        }

    _get_params_combine_gvcfs = _get_params_discover
    _get_params_genotype = _get_params_discover


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
        if self.report_per_donor and self.parent.effective_group_by == "cohort":
            for members in self.parent.cohort_members.values():
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

    published_keys = ("txt", "txt_md5")

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


class BcftoolsRohStepPart(ReportStepPart):
    """ROH calling with `bcftools roh`."""

    name = "bcftools_roh"
    report_per_donor = False
    published_keys = ("txt", "txt_md5")
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
        return self.parent.get_cohort_libraries(index_library_name)

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

    published_keys = ("bw", "bw_md5")

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
            "ignore_chroms": sorted(set(self.config.ignore_chroms)),
            "padding": self.config.mutect2.padding,
            "java_options": self.config.mutect2.contamination.java_options,
            "extra_arguments": self.config.mutect2.contamination.extra_arguments,
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
                    raise ValueError(
                        f"Normal sample for tumor {wildcards.library_name} required but not found."
                    )

        input_files.update(self.parent.reference_inputs())

        if self.config.depends_on.germline_resource:
            input_files["germline_resource"] = self.parent.get_upstream_paths(
                "germline_resource"
            ).vcf
        if self.config.depends_on.panel_of_normals:
            pon = self.parent.get_upstream_paths("panel_of_normals")
            input_files["panel_of_normals"] = pon.vcf
            input_files["panel_of_normals_tbi"] = pon.vcf_tbi

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
            "common_variants": self.parent.get_upstream_paths("common_variants").vcf,
        }

    def _get_input_files_pileup_tumor(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        return {
            "bam": alignments.bam,
            "bai": alignments.bai,
            **self.parent.reference_inputs(),
            "common_variants": self.parent.get_upstream_paths("common_variants").vcf,
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
        df = self.parent.build_library_dataframe()
        if df.empty:
            return

        selection = getattr(self.config, "library_selection", None)
        if selection:
            from snappy_pipeline.workflows.abstract import apply_library_selection

            df = apply_library_selection(df, selection, "cancer")
        else:
            df = df[df["role"] == "tumor"]

        tpl = "output/{library_name}/out/{library_name}.vcf.gz"
        for _, row in df.iterrows():
            for ext in ("", ".md5", ".tbi", ".tbi.md5"):
                yield tpl.format(library_name=row["library_name"]) + ext


class VariantCallingWorkflow(BaseStep):
    """Workflow implementation for germline variant calling"""

    name = "variant_calling"
    produces = [
        DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"})),
        DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel", "filtered"})),
    ]
    config_model_class = VariantCallingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """Germline small variants, or for mutect2 all somatic calls plus the PASS calls."""
        if config.tool != Tool.mutect2:
            return (DataSignature(DataType.VARIANTS, frozenset({"germline", "snv", "indel"})),)
        calls = DataSignature(DataType.VARIANTS, frozenset({"somatic", "snv", "indel"}))
        return (calls, calls.with_tags("filtered"))

    sheet_shortcut_class = GermlineCaseSheet

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
        if (
            self.config.depends_on.kits
            and self.config.bcftools_roh
            and self.config.bcftools_roh.enabled
        ):
            for index_library_name in self.output_entities:
                self.common_kit(
                    "kits", self.sub_steps["bcftools_roh"].vcf_libraries(index_library_name)
                )

    def get_result_files(self) -> list[str]:
        """The files of the caller and of the enabled reports"""
        names = [self.config.tool]
        for name in ("baf_file_generation", "bcftools_stats", "bcftools_roh"):
            name_config = self.config.get(name)
            if name_config and name_config.enabled:
                names.append(name)
        return [path for name in names for path in self.sub_steps[name].get_result_files()]

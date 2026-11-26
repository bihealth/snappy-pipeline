from typing import Any

from snakemake.io import directory
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import RuleOptions, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.common.sv_calling import VCF_OUTPUTS, SvCallingGetResultFilesMixin

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class MeltStepPart(SvCallingGetResultFilesMixin, BaseStepPart):
    """MEI calling using MELT

    We implement the workflow as per-pedigree calling.  Generally, this leads to consistent
    positions within each pedigree but not necessarily across the whole cohort.

    Note that MELT is not free software, so further setup is needed.
    """

    name = "melt"

    actions = (
        "preprocess",
        "indiv_analysis",
        "group_analysis",
        "genotype",
        "make_vcf",
        "merge_vcf",
    )

    wrappers = {
        "merge_vcf": "melt/merge_vcf",
        "make_vcf": "melt/make_vcf",
        "genotype": "melt/genotype",
        "group_analysis": "melt/group_analysis",
        "indiv_analysis": "melt/indiv_analysis",
        "preprocess": "melt/preprocess",
    }

    rule_options = {
        "preprocess": RuleOptions(constraints={"library_name": r"[^.]+", "me_type": r"[^.]+"}),
        "indiv_analysis": RuleOptions(constraints={"library_name": r"[^.]+", "me_type": r"[^.]+"}),
        "group_analysis": RuleOptions(constraints={"related_group": r"[^.]+", "me_type": r"[^.]+"}),
        "genotype": RuleOptions(constraints={"related_group": r"[^.]+", "me_type": r"[^.]+"}),
        "make_vcf": RuleOptions(constraints={"related_group": r"[^.]+"}),
        "merge_vcf": RuleOptions(constraints={"library_name": r"[^.]+"}),
    }

    #: The MELT commands of the other actions take no ``extra_arguments``
    extra_arguments_action = "indiv_analysis"

    intermediate_outputs = {
        "preprocess": ("orig_bam", "orig_bai", "disc_bam", "disc_bai"),
        "indiv_analysis": ("out",),
        "group_analysis": ("pre_geno",),
        "genotype": ("tsv",),
        "make_vcf": ("vcf", "vcf_tbi"),
    }

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="16GB")

    #: MELT calls only the pedigrees of the foreground libraries
    pedigrees_include_background = False

    @dictify
    def _get_input_files_preprocess(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar").items()

    @dictify
    def _get_output_files_preprocess(self):
        # MELT infers the sample name from the BAM file name instead of the BAM header.
        prefix = "work/{library_name}/out/{library_name}"
        yield "orig_bam", f"{prefix}.bam"
        yield "orig_bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield "disc_fq", f"{prefix}.bam.fq"

    @dictify
    def _get_log_file_preprocess(self):
        yield from self._get_log_file_with_prefix(
            "work/{library_name}/log/{library_name}.preprocess"
        ).items()

    @dictify
    def _get_input_files_indiv_analysis(self, wildcards):
        prefix = f"work/{wildcards.library_name}/out/{wildcards.library_name}"
        yield "orig_bam", f"{prefix}.bam"
        yield "orig_bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar", "me_refs").items()

    @dictify
    def _get_output_files_indiv_analysis(self):
        # IndivAnalysis names its files itself
        yield "out", directory("work/indiv_analysis.{library_name}.{me_type}/out")

    @dictify
    def _get_log_file_indiv_analysis(self):
        yield from self._get_log_file_with_prefix(
            "work/indiv_analysis.{library_name}.{me_type}/log/{library_name}.{me_type}"
        ).items()

    def _get_params_indiv_analysis(self, wildcards: Wildcards) -> dict[str, Any]:
        """IndivAnalysis takes MELT's ``-r``, ``-b`` and, for targeted data, ``-exome``"""
        params = self._get_params_group_analysis(wildcards)
        params["exclude_contigs"] = self.config.melt.exclude_contigs
        if self.parent.targeted:
            params["exome"] = True
        return params

    @dictify
    def _get_input_files_group_analysis(self, wildcards):
        yield from self.parent.reference_inputs().items()
        indiv_analysis = []
        for library_name in self.parent.pedigrees(include_background=True)[wildcards.related_group]:
            indiv_analysis.append(f"work/indiv_analysis.{library_name}.{wildcards.me_type}/out")
        yield "indiv_analysis", indiv_analysis
        yield from self._melt_files("jar", "me_refs", "genes").items()

    @dictify
    def _get_output_files_group_analysis(self):
        infix = "group_analysis.{related_group}.{me_type}"
        # Genotype and MakeVCF read the sites to genotype
        yield "pre_geno", f"work/{infix}/out/{{me_type}}.pre_geno.tsv"
        exts = (
            "bed.list",
            "hum.list",
            "master.bed",
            "merged.hum_breaks.sorted.bam",
            "merged.hum_breaks.sorted.bam.bai",
        )
        yield "_more", [f"work/{infix}/out/{{me_type}}.{ext}" for ext in exts]

    @dictify
    def _get_log_file_group_analysis(self):
        yield from self._get_log_file_with_prefix(
            "work/group_analysis.{related_group}.{me_type}/log/{related_group}.{me_type}"
        ).items()

    def _get_params_group_analysis(self, wildcards: Wildcards) -> dict[str, Any]:
        """GroupAnalysis takes MELT's ``-r``"""
        params = self._get_params_genotype(wildcards)
        params["read_length"] = self.config.melt.read_length
        return params

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield "pre_geno", self._pre_geno(wildcards)
        prefix = f"work/{wildcards.library_name}/out/{wildcards.library_name}"
        yield "bam", f"{prefix}.bam"
        yield "bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar", "me_refs").items()

    @dictify
    def _get_output_files_genotype(self):
        infix = "genotype.{related_group}.{me_type}"
        yield "tsv", f"work/{infix}/out/{{library_name}}.{{me_type}}.tsv"

    @dictify
    def _get_log_file_genotype(self):
        yield from self._get_log_file_with_prefix(
            "work/genotype.{related_group}.{me_type}/log/{library_name}.{me_type}"
        ).items()

    def _get_params_genotype(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "me_refs_infix": self.config.melt.me_refs_infix,
            "me_type": wildcards.me_type,
        }

    @dictify
    def _get_input_files_make_vcf(self, wildcards):
        yield "pre_geno", self._pre_geno(wildcards)
        infix = f"genotype.{wildcards.related_group}.{wildcards.me_type}"
        paths = [
            f"work/{infix}/out/{library_name}.{wildcards.me_type}.tsv"
            for library_name in self.parent.pedigrees(include_background=True)[
                wildcards.related_group
            ]
        ]
        yield "genotype", paths
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar", "me_refs").items()

    @dictify
    def _get_output_files_make_vcf(self):
        out_dir = "work/make_vcf.{related_group}.{me_type}/out"
        yield "list_txt", f"{out_dir}/list.txt"
        yield "vcf", f"{out_dir}/{{related_group}}.{{me_type}}.final_comp.vcf.gz"
        yield "vcf_tbi", f"{out_dir}/{{related_group}}.{{me_type}}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_log_file_make_vcf(self):
        yield from self._get_log_file_with_prefix(
            "work/make_vcf.{related_group}.{me_type}/log/{related_group}.{me_type}"
        ).items()

    def _get_params_make_vcf(self, wildcards: Wildcards) -> dict[str, Any]:
        """MakeVCF takes MELT's ``-j``"""
        params = self._get_params_genotype(wildcards)
        params["max_no_call_percent"] = self.config.melt.max_no_call_percent
        return params

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        vcfs = []
        for me_type in self.config.melt.me_types:
            out_dir = f"work/make_vcf.{wildcards.library_name}.{me_type}/out"
            vcfs.append(f"{out_dir}/{wildcards.library_name}.{me_type}.final_comp.vcf.gz")
        yield "vcf", vcfs
        yield "vcf_tbi", [f"{vcf}.tbi" for vcf in vcfs]

    def _get_output_files_merge_vcf(self):
        return dict(VCF_OUTPUTS)

    @dictify
    def _get_log_file_merge_vcf(self):
        yield from self._get_log_file_with_prefix(
            "work/{library_name}/log/{library_name}.merge_vcf"
        ).items()

    @dictify
    def _get_log_file_with_prefix(self, prefix: str):
        """Return dict of log files whose paths start with ``prefix``"""
        yield "log", f"{prefix}.log"

    def _melt_files(self, *keys: str) -> dict[str, str]:
        """Return the MELT jar, support files directory and genes file (``keys`` selects)."""
        files = {
            "jar": self.config.melt.jar_file,
            "me_refs": self.config.melt.me_refs_path,
            "genes": self.parent.get_upstream_paths("gene_regions").bed,
        }
        return {key: files[key] for key in keys}

    def _pre_geno(self, wildcards) -> str:
        """The sites of GroupAnalysis, which Genotype and MakeVCF find in its directory (``-p``)"""
        infix = f"group_analysis.{wildcards.related_group}.{wildcards.me_type}"
        return f"work/{infix}/out/{wildcards.me_type}.pre_geno.tsv"

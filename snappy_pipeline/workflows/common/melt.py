from typing import Any

from snakemake.io import touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.common.sv_calling import SvCallingGetResultFilesMixin

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class MeltStepPart(SvCallingGetResultFilesMixin, BaseStepPart):
    """MEI calling using MELT

    We implement the workflow as per-pedigree calling.  Generally, this leads to consistent
    positions within each pedigree but not necessarily across the whole cohort.

    Note that MELT is not free software, so further setup is needed.
    """

    name = "melt"

    #: MELT calls only the pedigrees of foreground data sets
    pedigrees_include_background = False

    actions = (
        "preprocess",
        "indiv_analysis",
        "group_analysis",
        "genotype",
        "make_vcf",
        "merge_vcf",
    )

    intermediate_outputs = {
        "preprocess": ("orig_bam", "orig_bai", "disc_bam", "disc_bai"),
        "indiv_analysis": ("done",),
        "group_analysis": ("done",),
        "genotype": ("done",),
        "make_vcf": ("vcf", "vcf_tbi"),
    }

    #: The MELT commands of the other actions take no ``extra_arguments``
    extra_arguments_action = "indiv_analysis"

    default_resource_usage = ResourceUsage(threads=1, runtime="1d", mem="16GB")

    def _melt_files(self, *keys: str) -> dict[str, str]:
        """Return the MELT jar, support files directory and genes file (``keys`` selects)."""
        files = {
            "jar": self.config.melt.jar_file,
            "me_refs": self.config.melt.me_refs_path,
            "genes": self.parent.get_upstream_paths("gene_regions").bed,
        }
        return {key: files[key] for key in keys}

    @dictify
    def _get_log_file_with_prefix(self, prefix: str):
        """Return dict of log files whose paths start with ``prefix``"""
        yield "log", f"{prefix}.log"

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
        yield "done", touch("work/indiv_analysis.{library_name}.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_log_file_indiv_analysis(self):
        yield from self._get_log_file_with_prefix(
            "work/indiv_analysis.{library_name}.{me_type}/log/{library_name}.{me_type}"
        ).items()

    @dictify
    def _get_input_files_group_analysis(self, wildcards):
        yield from self.parent.reference_inputs().items()
        indiv_analysis = []
        for library_name in self.parent.pedigrees()[wildcards.index_library_name]:
            infix = f"indiv_analysis.{library_name}.{wildcards.me_type}"
            indiv_analysis.append(f"work/{infix}/out/.done.{library_name}")
        yield "indiv_analysis", indiv_analysis
        yield from self._melt_files("jar", "me_refs", "genes").items()

    @dictify
    def _get_output_files_group_analysis(self):
        infix = "group_analysis.{index_library_name}.{me_type}"
        yield "done", touch(f"work/{infix}/out/.done")
        exts = (
            "bed.list",
            "hum.list",
            "master.bed",
            "merged.hum_breaks.sorted.bam",
            "merged.hum_breaks.sorted.bam.bai",
            "pre_geno.tsv",
        )
        yield "_more", [f"work/{infix}/out/{{me_type}}.{ext}" for ext in exts]

    @dictify
    def _get_log_file_group_analysis(self):
        yield from self._get_log_file_with_prefix(
            "work/group_analysis.{index_library_name}.{me_type}/log/{index_library_name}.{me_type}"
        ).items()

    @dictify
    def _get_input_files_genotype(self, wildcards):
        infix_done = f"group_analysis.{wildcards.index_library_name}.{wildcards.me_type}"
        yield "done", f"work/{infix_done}/out/.done"
        prefix = f"work/{wildcards.library_name}/out/{wildcards.library_name}"
        yield "bam", f"{prefix}.bam"
        yield "bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar", "me_refs").items()

    @dictify
    def _get_output_files_genotype(self):
        infix = "genotype.{index_library_name}.{me_type}"
        yield "done", touch(f"work/{infix}/out/.done.{{library_name}}")
        yield "_more", [f"work/{infix}/out/{{library_name}}.{{me_type}}.tsv"]

    @dictify
    def _get_log_file_genotype(self):
        yield from self._get_log_file_with_prefix(
            "work/genotype.{index_library_name}.{me_type}/log/{library_name}.{me_type}"
        ).items()

    @dictify
    def _get_input_files_make_vcf(self, wildcards):
        infix = f"group_analysis.{wildcards.index_library_name}.{wildcards.me_type}"
        yield "group_analysis", f"work/{infix}/out/.done"
        infix = f"genotype.{wildcards.index_library_name}.{wildcards.me_type}"
        paths = [
            f"work/{infix}/out/.done.{library_name}"
            for library_name in self.parent.pedigrees()[wildcards.index_library_name]
        ]
        yield "genotype", paths
        yield from self.parent.reference_inputs().items()
        yield from self._melt_files("jar", "me_refs").items()

    @dictify
    def _get_log_file_make_vcf(self):
        yield from self._get_log_file_with_prefix(
            "work/make_vcf.{index_library_name}.{me_type}/log/{index_library_name}.{me_type}"
        ).items()

    @dictify
    def _get_output_files_make_vcf(self):
        out_dir = "work/make_vcf.{index_library_name}.{me_type}/out"
        yield "list_txt", f"{out_dir}/list.txt"
        yield "done", touch(f"{out_dir}/.done")
        yield "vcf", f"{out_dir}/{{index_library_name}}.{{me_type}}.final_comp.vcf.gz"
        yield "vcf_tbi", f"{out_dir}/{{index_library_name}}.{{me_type}}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        vcfs = []
        for me_type in self.config.melt.me_types:
            out_dir = f"work/make_vcf.{wildcards.library_name}.{me_type}/out"
            vcfs.append(f"{out_dir}/{wildcards.library_name}.{me_type}.final_comp.vcf.gz")
        yield "vcf", vcfs
        yield "vcf_tbi", [f"{vcf}.tbi" for vcf in vcfs]

    @dictify
    def _get_output_files_merge_vcf(self):
        prefix = "work/{library_name}/out/{library_name}"
        work_files = {
            "vcf": f"{prefix}.vcf.gz",
            "vcf_md5": f"{prefix}.vcf.gz.md5",
            "vcf_tbi": f"{prefix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"{prefix}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()

    @dictify
    def _get_log_file_merge_vcf(self):
        yield from self._get_log_file_with_prefix(
            "work/{library_name}/log/{library_name}.merge_vcf"
        ).items()

    def _get_params_group_analysis(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "me_refs_infix": self.config.melt.me_refs_infix,
            "me_type": wildcards.me_type,
        }

    _get_params_genotype = _get_params_group_analysis
    _get_params_make_vcf = _get_params_group_analysis

    def _get_params_indiv_analysis(self, wildcards: Wildcards) -> dict[str, Any]:
        """Only IndivAnalysis takes MELT's ``-exome``"""
        params = self._get_params_group_analysis(wildcards)
        if self.parent.name == "sv_calling_targeted":
            params["exome"] = True
        return params

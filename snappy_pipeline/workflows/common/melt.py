import re
from itertools import chain
from typing import Any

from biomedsheets.shortcuts import is_not_background
from snakemake.io import touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.common import (
    ForwardResourceUsageMixin,
    ForwardSnakemakeFilesMixin,
)
from snappy_pipeline.workflows.common.sv_calling import SvCallingGetResultFilesMixin

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


class MeltStepPart(
    SvCallingGetResultFilesMixin,
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    BaseStepPart,
):
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

    _resource_usage = ResourceUsage(
        threads=1,
        runtime="1d",
        mem="16GB",
    )
    resource_usage_dict = {
        "preprocess": _resource_usage,
        "indiv_analysis": _resource_usage,
        "group_analysis": _resource_usage,
        "genotype": _resource_usage,
        "make_vcf": _resource_usage,
        "merge_vcf": _resource_usage,
        "reorder_vcf": _resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)
        #: All individual's primary NGS libraries
        self.all_dna_ngs_libraries = []
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    self.all_dna_ngs_libraries.append(donor.dna_ngs_library.name)
        #: Linking NGS libraries to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def _get_log_file_with_prefix(self, prefix: str):
        """Return dict of log files whose paths start with ``prefix``"""
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield key + "_md5", f"{prefix}{ext}.md5"

    @dictify
    def _get_input_files_preprocess(self, wildcards):
        alignments = self.parent.get_upstream_paths(
            "alignments", library_name=wildcards.library_name
        )
        yield "bam", alignments.bam
        yield "bai", alignments.bai
        yield "reference", self.w_config.static_data_config.reference.path

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
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "reference", self.w_config.static_data_config.reference.path

    @dictify
    def _get_output_files_indiv_analysis(self):
        yield "done", touch("work/indiv_analysis.{library_name}.{me_type}/out/.done.{library_name}")

    @dictify
    def _get_log_file_indiv_analysis(self):
        yield from self._get_log_file_with_prefix(
            "work/indiv_analysis.{library_name}.{me_type}/log/{library_name}.{me_type}"
        ).items()

    @listify
    def _get_input_files_group_analysis(self, wildcards):
        yield self.w_config.static_data_config.reference.path
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        for member in pedigree.donors:
            if member.dna_ngs_library:
                library_name = member.dna_ngs_library.name
                infix = f"indiv_analysis.{library_name}.{wildcards.me_type}"
                yield f"work/{infix}/out/.done.{library_name}"

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
        yield "bam", f"work/{wildcards.library_name}/out/{wildcards.library_name}.bam"
        yield "reference", self.w_config.static_data_config.reference.path

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
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        paths = []
        for member in pedigree.donors:
            if member.dna_ngs_library:
                infix = f"genotype.{wildcards.index_library_name}.{wildcards.me_type}"
                paths.append(f"work/{infix}/out/.done.{member.dna_ngs_library.name}")
        yield "genotype", paths
        yield "reference", self.w_config.static_data_config.reference.path

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
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(work_files.values(), self.get_log_file("merge_vcf").values())
            ],
        )

    @dictify
    def _get_log_file_merge_vcf(self):
        yield from self._get_log_file_with_prefix(
            "work/{library_name}/log/{library_name}.merge_vcf"
        ).items()

    def _get_params_preprocess(self, wildcards: Wildcards) -> dict[str, Any]:
        params = {
            "config": self.config.melt.model_dump(by_alias=True),
        }
        if self.parent.name == "sv_calling_targeted":
            params["exome"] = True
        if getattr(wildcards, "me_type", None):
            params["me_type"] = getattr(wildcards, "me_type")
        return params

    _get_params_indiv_analysis = _get_params_preprocess
    _get_params_group_analysis = _get_params_preprocess
    _get_params_genotype = _get_params_preprocess
    _get_params_make_vcf = _get_params_preprocess
    _get_params_merge_vcf = _get_params_preprocess

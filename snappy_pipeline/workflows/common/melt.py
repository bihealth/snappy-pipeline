from itertools import chain
import re
import typing

from snakemake.io import touch

from biomedsheets.shortcuts import is_not_background
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
        time="1-00:00:00",
        memory="16G",
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
    def _get_log_file_with_infix(self, infix: str, *, suffix: typing.Optional[str] = None):
        """Return dict of log files in the "log" directory"""
        if suffix:
            suffix_str = f"_{suffix}"
        else:
            suffix_str = ""
        prefix = f"work/{infix}/log/{infix}.sv_calling{suffix_str}"
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
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")
        yield "bai", ngs_mapping(f"output/{infix}/out/{infix}.bam.bai")

    @dictify
    def _get_output_files_preprocess(self):
        # Note that mapper is not part of the output BAM file as MELT infers sample file from BAM
        # file name instead of using sample name from BAM header.
        prefix = "work/{mapper}.melt_preprocess.{library_name}/out/{library_name}"
        yield "orig_bam", f"{prefix}.bam"
        yield "orig_bai", f"{prefix}.bam.bai"
        yield "disc_bam", f"{prefix}.bam.disc"
        yield "disc_bai", f"{prefix}.bam.disc.bai"
        yield "disc_fq", f"{prefix}.bam.fq"

    @dictify
    def _get_log_file_preprocess(self):
        yield from self._get_log_file_with_infix("{mapper}.melt_preprocess.{library_name}").items()

    @dictify
    def _get_input_files_indiv_analysis(self, wildcards):
        infix = f"{wildcards.mapper}.melt_preprocess.{wildcards.library_name}"
        yield "orig_bam", f"work/{infix}/out/{wildcards.library_name}.bam"
        yield "disc_bam", f"work/{infix}/out/{wildcards.library_name}.bam.disc"

    @dictify
    def _get_output_files_indiv_analysis(self):
        infix = "{mapper}.melt_indiv_analysis.{library_name}.{me_type}"
        yield "done", touch(f"work/{infix}/out/.done.{{library_name}}")

    @dictify
    def _get_log_file_indiv_analysis(self):
        yield from self._get_log_file_with_infix(
            "{mapper}.melt_indiv_analysis.{library_name}.{me_type}", suffix="_{library_name}"
        ).items()

    @listify
    def _get_input_files_group_analysis(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        for member in pedigree.donors:
            if member.dna_ngs_library:
                infix = f"{wildcards.mapper}.melt_indiv_analysis.{member.dna_ngs_library.name}.{wildcards.me_type}"
                yield f"work/{infix}/out/.done.{member.dna_ngs_library.name}"

    @dictify
    def _get_output_files_group_analysis(self):
        infix = "{mapper}.melt_group_analysis.{index_library_name}.{me_type}"
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
        yield from self._get_log_file_with_infix(
            "{mapper}.melt_group_analysis.{index_library_name}.{me_type}"
        ).items()

    @dictify
    def _get_input_files_genotype(self, wildcards):
        infix_done = f"{wildcards.mapper}.melt_group_analysis.{wildcards.index_library_name}.{wildcards.me_type}"
        yield "done", f"work/{infix_done}/out/.done".format(**wildcards)
        infix_bam = f"{wildcards.mapper}.melt_preprocess.{wildcards.library_name}"
        yield "bam", f"work/{infix_bam}/out/{wildcards.library_name}.bam"

    @dictify
    def _get_output_files_genotype(self):
        infix = "{mapper}.melt_genotype.{index_library_name}.{me_type}"
        yield "done", touch(f"work/{infix}/out/.done.{{library_name}}")
        yield "_more", [f"work/{infix}/out/{{library_name}}.{{me_type}}.tsv"]

    @dictify
    def _get_log_file_genotype(self):
        yield from self._get_log_file_with_infix(
            "{mapper}.melt_genotype.{index_library_name}.{me_type}", suffix="_{library_name}"
        ).items()

    @dictify
    def _get_input_files_make_vcf(self, wildcards):
        infix = f"{wildcards.mapper}.melt_group_analysis.{wildcards.index_library_name}.{wildcards.me_type}"
        yield "group_analysis", f"work/{infix}/out/.done"
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_library_name]
        paths = []
        for member in pedigree.donors:
            if member.dna_ngs_library:
                infix = f"{wildcards.mapper}.melt_genotype.{wildcards.index_library_name}.{wildcards.me_type}"
                paths.append(f"work/{infix}/out/.done.{member.dna_ngs_library.name}")
        yield "genotype", paths

    @dictify
    def _get_log_file_make_vcf(self):
        yield from self._get_log_file_with_infix(
            "{mapper}.melt_make_vcf.{index_library_name}.{me_type}"
        ).items()

    @dictify
    def _get_output_files_make_vcf(self):
        infix = "{mapper}.melt_make_vcf.{index_library_name}.{me_type}"
        yield "list_txt", f"work/{infix}/out/list.txt"
        yield "done", touch(f"work/{infix}/out/.done")
        yield "vcf", f"work/{infix}/out/{infix}.final_comp.vcf.gz"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.final_comp.vcf.gz.tbi"

    @dictify
    def _get_input_files_merge_vcf(self, wildcards):
        vcfs = []
        for me_type in self.config.melt.me_types:
            infix = f"{wildcards.mapper}.melt_make_vcf.{wildcards.library_name}.{me_type}"
            vcfs.append(f"work/{infix}/out/{infix}.final_comp.vcf.gz")
        yield "vcf", vcfs

    @dictify
    def _get_output_files_merge_vcf(self):
        infix = "{mapper}.melt.{library_name}"
        work_files = {
            "vcf": f"work/{infix}/out/{infix}.vcf.gz",
            "vcf_md5": f"work/{infix}/out/{infix}.vcf.gz.md5",
            "vcf_tbi": f"work/{infix}/out/{infix}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{infix}/out/{infix}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("merge_vcf").values())
        ]

    @dictify
    def _get_log_file_merge_vcf(self):
        yield from self._get_log_file_with_infix("{mapper}.melt.{library_name}").items()

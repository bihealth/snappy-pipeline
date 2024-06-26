"""Implementation of the ``sv_calling_wgs`` step"""

import re
from itertools import chain

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import (
    ForwardResourceUsageMixin,
    ForwardSnakemakeFilesMixin,
)
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart
from snappy_pipeline.workflows.common.sv_calling import (
    SvCallingGetLogFileMixin,
    SvCallingGetResultFilesMixin,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_wrappers.tools.genome_windows import yield_regions

from .model import SvCallingWgs as SvCallingWgsConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Available (short) DNA WGS SV callers
DNA_WGS_SV_CALLERS = ("delly2", "manta", "popdel", "melt", "gcnv")

#: Available (long) DNA WGS SV callers
LONG_DNA_WGS_SV_CALLERS = ("pb_honey_spots", "sniffles", "sniffles2")

#: Default configuration for the sv_calling_wgs step
DEFAULT_CONFIG = SvCallingWgsConfigModel.default_config_yaml_string()


class GcnvWgsStepPart(RunGcnvStepPart):
    """WGS CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        # No mapping given for WGS; we will use the "wgs" one for all.
        for donor in self.parent.all_donors():
            if donor.dna_ngs_library:
                yield donor.dna_ngs_library.name, "wgs"


def escape_dots_dashes(s: str) -> str:
    """Escape dots and dashes with double-underscore constructs."""
    return s.replace("_", "__under__").replace(".", "__dot__").replace("-", "__hyphen__")


class PopDelStepPart(
    SvCallingGetResultFilesMixin,
    SvCallingGetLogFileMixin,
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    BaseStepPart,
):
    """WGS SV identification using PopDel.

    PopDel requires whole-cohort wide calling.  We implement cohort-wise calling to speed
    things up a bit more.
    """

    name = "popdel"
    actions = ("profile", "call", "concat_calls", "reorder_vcf")

    _resource_usage = ResourceUsage(
        threads=2,
        time="4-00:00:00",
        memory="24G",
    )
    resource_usage_dict = {
        "profile": _resource_usage,
        "call": _resource_usage,
        "concat_calls": _resource_usage,
        "reorder_vcf": _resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Build shortcut from library name to library info
        self.library_name_to_library = {}
        for sheet in self.parent.shortcut_sheets:
            self.library_name_to_library.update(sheet.library_name_to_library)

    def get_library_extra_infos(self, wildcards):
        """Returns library extra infos for the given library name"""
        return self.library_name_to_library[wildcards.library_name].ngs_library.extra_infos

    @dictify
    def _get_input_files_profile(self, wildcards):
        """Return input files for "call" action"""
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def _get_output_files_profile(self):
        infix = "{mapper}.popdel_profile.{library_name}"
        yield "profile", f"work/{infix}/out/{infix}.profile"
        yield "profile_md5", f"work/{infix}/out/{infix}.profile.md5"

    @dictify
    def _get_input_files_call(self, wildcards):
        paths = []
        for donor in self._donors_with_dna_ngs_library():
            infix = f"{wildcards.mapper}.popdel_profile.{donor.dna_ngs_library.name}"
            paths.append(f"work/{infix}/out/{infix}.profile")
        yield "profile", paths

    def _donors_with_dna_ngs_library(self):
        """Yield donors with DNA NGS library"""
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                if donor.dna_ngs_library:
                    yield donor

    @dictify
    def _get_output_files_call(self):
        infix = self._get_log_file_infix_call()
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    def _get_log_file_infix_call(self):
        return "{mapper}.popdel_call.{chrom}-{begin}-{end}"

    @dictify
    def _get_input_files_concat_calls(self, wildcards):
        window_size = self.config.popdel.window_size
        padding = self.config.popdel.max_sv_size
        vcfs = []
        with open(self._get_fai_path(), "rt") as fai_file:
            for r in yield_regions(
                fai_file,
                window_size=window_size,
                padding=padding,
                ignore_chroms=self._get_ignore_chroms(),
            ):
                if r.begin == 0:
                    r.begin = 1
                chrom = escape_dots_dashes(r.chrom)
                infix = f"{wildcards.mapper}.popdel_call.{chrom}-{r.begin}-{r.end}"
                vcfs.append(f"work/{infix}/out/{infix}.vcf.gz")
        yield "vcf", vcfs

    def _get_fai_path(self):
        return self.w_config.static_data_config.reference.path + ".fai"

    def _get_ignore_chroms(self):
        return self.config.ignore_chroms

    @dictify
    def _get_output_files_concat_calls(self):
        infix = self._get_log_file_infix_concat_calls()
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"
        yield "vcf_md5", f"work/{infix}/out/{infix}.vcf.gz.md5"
        yield "vcf_tbi", f"work/{infix}/out/{infix}.vcf.gz.tbi"
        yield "vcf_tbi_md5", f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"

    def _get_log_file_infix_concat_calls(self):
        return "{mapper}.popdel_concat_calls"

    @dictify
    def _get_input_files_reorder_vcf(self, wildcards):
        infix = f"{wildcards.mapper}.popdel_concat_calls"
        yield "vcf", f"work/{infix}/out/{infix}.vcf.gz"

    @dictify
    def _get_output_files_reorder_vcf(self):
        infix = "{mapper}.popdel.{library_name}"
        work_files = {}
        work_files["vcf"] = f"work/{infix}/out/{infix}.vcf.gz"
        work_files["vcf_md5"] = f"work/{infix}/out/{infix}.vcf.gz.md5"
        work_files["vcf_tbi"] = f"work/{infix}/out/{infix}.vcf.gz.tbi"
        work_files["vcf_tbi_md5"] = f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"
        yield from work_files.items()
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(
                    work_files.values(), self.get_log_file("reorder_vcf").values()
                )
            ],
        )

    def get_ped_members(self, wildcards):
        """Used in Snakefile to rule ``sv_calling_wgs_popdel_reorder_vcf``"""
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        return " ".join(
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        )


class Sniffles2StepPart(BaseStepPart):
    """WGS SV identification using Sniffles 2"""

    name = "sniffles2"
    actions = ("bam_to_snf", "snf_to_vcf")

    _resource_usage = ResourceUsage(threads=2, time="0-02:00:00", memory="4G")
    resource_usage_dict = {
        "bam_to_snf": _resource_usage,
        "snf_to_vcf": _resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = (
            "work/{mapper}.sniffles2.{index_ngs_library}/out/"
            "{mapper}.sniffles2.{index_ngs_library}{ext}"
        )
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_bam_to_snf(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        infix = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{infix}/out/{infix}.bam")

    @dictify
    def _get_output_files_bam_to_snf(self):
        infix = "{mapper}.sniffles2_bam_to_snf.{library_name}"
        yield "snf", f"work/{infix}/out/{infix}.snf"

    @dictify
    def _get_input_files_snf_to_vcf(self, wildcards):
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        snfs = []
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.sniffles2_bam_to_snf.{donor.dna_ngs_library.name}.snf"
                snfs.append(f"work/{infix}/out/{infix}.snf")
        yield "snf", snfs

    @dictify
    def _get_output_files_snf_to_vcf(self):
        infix = "{mapper}.sniffles2.{index_ngs_library}"
        yield "snf", f"work/{infix}/out/{infix}.snf"


class SvCallingWgsWorkflow(BaseStep):
    """Perform (germline) WGS SV calling"""

    name = "sv_calling_wgs"
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SvCallingWgsConfigModel,
            previous_steps=(NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                Delly2StepPart,
                MantaStepPart,
                PopDelStepPart,
                GcnvWgsStepPart,
                MeltStepPart,
                # Sniffles2StepPart,
                WritePedigreeStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

    @listify
    def all_donors(self, include_background=True):
        """Return list of all donors in sample sheet."""
        sheets = self.shortcut_sheets
        if not include_background:
            filter(is_not_background, sheets)
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample sheets.
        """
        for sub_step in self.sub_steps.values():
            yield from sub_step.get_result_files()

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )

import random
from typing import Any

from snakemake.io import Wildcards
from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import ResourceUsage
from snappy_pipeline.workflows.common.samplesheet import tumor_to_normal_mapping

from snappy_pipeline.workflows.any_variant_filtration import (
    AnyVariantFiltrationWorkflow,
    OneFilterWithBamStepPart,
)

from snappy_pipeline.workflows.any_variant_calling.model import VariantOrigin
from .model import SomaticVariantFiltration as SomaticVariantFiltrationConfigModel
from .model import Ebfilter as EbfilterConfig


class OneFilterEbfilterStepPart(OneFilterWithBamStepPart):
    name = "one_ebfilter"
    filter_name = "ebfilter"

    #: Class available actions
    actions = ("run", "write_panel")

    resource_usage = {
        "run": ResourceUsage(threads=1, time="24:00:00", memory=f"{2 * 1024}M"),
        "write_panel": ResourceUsage(threads=1, time="01:00:00", memory=f"{2 * 1024}M"),
    }

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path to input or previous filter vcf file & normal/tumor bams"""
        parent = super(OneFilterEbfilterStepPart, self)._get_input_files_run
        yield from parent(wildcards).items()
        cfg: EbfilterConfig = self._get_args(wildcards)
        sample_files = cfg["path_panel_of_normals_sample_list"]
        if not sample_files:
            sample_files = self._get_output_files_write_panel()["txt"].format(**wildcards)
        yield "txt", sample_files

    def _get_output_files_write_panel(self):
        return {
            "txt": "work/{mapper}.eb_filter.panel_of_normals/out/{mapper}.eb_filter.panel_of_normals.txt"
        }

    def get_output_files(self, action):
        output_files = super(OneFilterEbfilterStepPart, self).get_output_files(action)
        if action == "write_panel":
            output_files = self._get_output_files_write_panel()
        return output_files

    def _get_args(self, wildcards: Wildcards) -> dict[str, Any]:
        """Return dkfz parameters to parameters"""
        return super(OneFilterEbfilterStepPart, self)._get_args(wildcards) | {
            "has_annotation": self.config.has_annotation,
        }

    def write_panel_of_normals_file(self, wildcards):
        """Write out file with paths to panels-of-normal"""
        output_path = self.get_output_files("write_panel")["txt"].format(**wildcards)
        with open(output_path, "wt") as outf:
            for bam_path in self._get_panel_of_normal_bams(wildcards):
                print(bam_path, file=outf)

    @listify
    def _get_panel_of_normal_bams(self, wildcards):
        """Return list of "panel of normal" BAM files."""
        libraries = list(
            filter(
                lambda x: x is not None and x != "", self.parent.tumor_to_normal_mapping.values()
            )
        )
        libraries.sort()

        for filter_cfg in self.config.filter_list:
            if list(filter_cfg.keys())[0] == "ebfilter":
                cfg: EbfilterConfig = list(filter_cfg.values())[0]
                break
        random.seed(cfg.shuffle_seed)
        lib_count = cfg["panel_of_normals_size"]
        random.shuffle(libraries)
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}"
        for library in libraries[:lib_count]:
            yield ngs_mapping(tpl.format(normal_library=library, **wildcards) + ".bam")


class SomaticVariantFiltrationWorkflow(AnyVariantFiltrationWorkflow):
    name = "somatic_variant_filtration"
    variant_origin = VariantOrigin.SOMATIC
    model_class = SomaticVariantFiltrationConfigModel

    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=False, allow_missing_tumor=False)
    }

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(workflow, config, config_lookup_paths, config_paths, workdir)

        # Register sub step classes so the sub steps are available
        sub_steps = list(map(lambda x: x.__class__, self.sub_steps.values()))
        self.register_sub_step_classes(sub_steps + [OneFilterEbfilterStepPart])

        self.tumor_to_normal_mapping = tumor_to_normal_mapping(self.table)

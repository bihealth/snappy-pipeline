import os
import sys
from collections import OrderedDict

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_annotation import (
    ANNOTATION_TOOLS,
    SomaticVariantAnnotationWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)
from snappy_pipeline.workflows.somatic_variant_filtration import SomaticVariantFiltrationWorkflow

from .model import TumorMutationalBurden as TumorMutationalBurdenConfigModel

#: Extensions of files to create as main payload
EXT_VALUES = (".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("json", "json_md5")

#: Default configuration for the tmb calculation step
DEFAULT_CONFIG = TumorMutationalBurdenConfigModel.default_config_yaml_string()


class TumorMutationalBurdenCalculationStepPart(BaseStepPart):
    """Calculation tumor mutational burden for each sample"""

    name = "tmb_gathering"

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.config = parent.w_config.step_config["tumor_mutational_burden"]
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            # update function of OrderedDict
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )
        # Build mapping from donor name to donor.
        self.donors = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                self.donors[donor.name] = donor

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)

        additional_steps = ""
        if self.config.has_annotation:
            additional_steps += ".{anno_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                additional_steps += ".filtered"
            else:
                additional_steps += ".dkfz_bias_filter.eb_filter"
        base_name = "{mapper}.{var_caller}" + additional_steps + ".{tumor_library}"
        if self.config.filters:
            base_name += ".{filter}"
        if self.config.filtered_regions:
            base_name += ".{region}"
        tpl = os.path.join("output", base_name, "out", base_name)

        key_ext = {"vcf": ".vcf.gz", "vcf_tbi": ".vcf.gz.tbi"}
        variant_path = self.parent.sub_workflows["somatic_variant"]
        for key, ext in key_ext.items():
            yield key, variant_path(tpl + ext)

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)

        additional_steps = ""
        if self.config.has_annotation:
            additional_steps += ".{anno_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                additional_steps += ".filtered"
            else:
                additional_steps += ".dkfz_bias_filter.eb_filter"
        base_name = "{mapper}.{var_caller}" + additional_steps + ".tmb.{tumor_library}"
        if self.config.filters:
            base_name += ".{filter}"
        if self.config.filtered_regions:
            base_name += ".{region}"
        tpl = os.path.join("output", base_name, "out", base_name)

        key_ext = {"json": ".json"}
        for key, ext in key_ext.items():
            yield key, tpl + ext
            yield key + "_md5", tpl + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)

        additional_steps = ""
        if self.config.has_annotation:
            additional_steps += ".{anno_caller}"
        if self.config.is_filtered:
            if len(self.config.filters) == 0:
                additional_steps += ".filtered"
            else:
                additional_steps += ".dkfz_bias_filter.eb_filter"
        base_name = "{mapper}.{var_caller}" + additional_steps + ".tmb.{tumor_library}"
        if self.config.filters:
            base_name += ".{filter}"
        if self.config.filtered_regions:
            base_name += ".{region}"
        tpl = os.path.join("output", base_name, "log", base_name)

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, tpl + ext

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        mem_mb = 4 * 1024  # 4GB
        return ResourceUsage(
            threads=2,
            time="1:00:00",  # 1 hour
            memory=f"{mem_mb}M",
        )

    def get_params(self, action):
        self._validate_action(action)
        return getattr(self, "_get_params_run")

    def _get_params_run(self, wildcards):
        return {"missense_re": self.w_config.step_config["tumor_mutational_burden"].missense_regex}


class TumorMutationalBurdenCalculationWorkflow(BaseStep):
    """Perform TMB calculation"""

    name = "tumor_mutational_burden"
    sheet_shortcut_class = CancerCaseSheet
    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one."""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=TumorMutationalBurdenConfigModel,
            previous_steps=(
                SomaticVariantCallingWorkflow,
                SomaticVariantAnnotationWorkflow,
                SomaticVariantFiltrationWorkflow,
                NgsMappingWorkflow,
            ),
        )
        # Register sub workflows
        config = self.config
        sub_workflow = "somatic_variant_calling"
        if config.has_annotation:
            sub_workflow = "somatic_variant_annotation"
        if config.is_filtered:
            sub_workflow = "somatic_variant_filtration"
        self.register_sub_workflow(sub_workflow, config.path_somatic_variant, "somatic_variant")
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not config.tools_ngs_mapping:
            config.tools_ngs_mapping = self.w_config.step_config["ngs_mapping"].tools.dna
        if not config.tools_somatic_variant_calling:
            config.tools_somatic_variant_calling = self.w_config.step_config[
                "somatic_variant_calling"
            ].tools
        if not config.tools_somatic_variant_annotation:
            config.tools_somatic_variant_annotation = self.w_config.step_config[
                "somatic_variant_annotation"
            ].tools
        if config.is_filtered:
            if len(self.w_config.step_config["somatic_variant_filtration"].filter_list) > 0:
                config.filters = []
                config.filtered_regions = []
            else:
                if not config.filters:
                    config.filters = list(
                        self.w_config.step_config["somatic_variant_filtration"].filter_sets.keys()
                    )
                    config.filters.append("no_filter")
                if not config.filtered_regions:
                    config.filtered_regions = list(
                        self.w_config.step_config["somatic_variant_filtration"].exon_lists.keys()
                    )
                    config.filtered_regions.append("genome_wide")
        # Register sub step classes so the sub steps are available
        self.w_config.step_config["tumor_mutational_burden"] = config
        self.register_sub_step_classes((TumorMutationalBurdenCalculationStepPart, LinkOutStepPart))

    @listify
    def get_result_files(self):
        config = self.w_config.step_config["tumor_mutational_burden"]
        name_pattern = "{mapper}.{caller}"
        if config.has_annotation:
            name_pattern += ".{anno_caller}"
        if config.is_filtered:
            if len(config.filters) > 0:
                name_pattern += ".dkfz_bias_filter.eb_filter"
            else:
                name_pattern += ".filtered"
        name_pattern += ".tmb.{tumor_library.name}"
        if config.is_filtered and len(config.filters) > 0:
            name_pattern += ".{filter}.{region}"

        mappers = set(config.tools_ngs_mapping) & set(
            self.w_config.step_config["ngs_mapping"].tools.dna
        )
        assert len(mappers) > 0, "No valid mapper"
        callers = set(config.tools_somatic_variant_calling) & set(SOMATIC_VARIANT_CALLERS_MATCHED)
        assert len(callers) > 0, "No valid somatic variant caller"
        if config.has_annotation:
            anno_callers = set(config.tools_somatic_variant_annotation) & set(ANNOTATION_TOOLS)
            assert len(anno_callers) > 0, "No valid somatic variant annotation tool"
        else:
            anno_callers = []
        if config.is_filtered:
            filters = list(
                self.w_config.step_config["somatic_variant_filtration"].filter_sets.keys()
            )
            filters.append("no_filter")
            filters = set(filters) & set(config.filters)
            regions = list(
                self.w_config.step_config["somatic_variant_filtration"].exon_lists.keys()
            )
            regions.append("genome_wide")
            regions = set(regions) & set(config.filtered_regions)
        else:
            filters = []
            regions = []

        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=mappers,
            caller=callers,
            anno_caller=anno_callers,
            filter=filters,
            region=regions,
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=mappers,
            caller=callers,
            anno_caller=anno_callers,
            filter=filters,
            region=regions,
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

    def _yield_result_files_matched(self, tpl, **kwargs):
        """Build output paths from path template and extension list.

        This function returns the results from the matched somatic variant callers such as
        Mutect.
        """
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for sample_pair in sheet.all_sample_pairs:
                if (
                    not sample_pair.tumor_sample.dna_ngs_library
                    or not sample_pair.normal_sample.dna_ngs_library
                ):
                    msg = (
                        "INFO: sample pair for cancer bio sample {} has is missing primary"
                        "normal or primary cancer NGS library"
                    )
                    print(msg.format(sample_pair.tumor_sample.name), file=sys.stderr)
                    continue
                yield from expand(
                    tpl,
                    tumor_library=[sample_pair.tumor_sample.dna_ngs_library],
                    **kwargs,
                )

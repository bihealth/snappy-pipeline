import os
import sys

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, LinkOutStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow, ResourceUsage
from snappy_pipeline.workflows.somatic_variant_calling import (
    SOMATIC_VARIANT_CALLERS_MATCHED,
    SomaticVariantCallingWorkflow,
)

#: Extensions of files to create as main payload
EXT_VALUES = (".json", ".json.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("json", "json_md5")

#: Default configuration for the tmb calculation step
DEFAULT_CONFIG = r"""
step_config:
    somatic_variant_checking:
        path_somatic_variant_calling: ../somatic_variant_calling   # REQUIRED
        tools_ngs_mapping: []      # default to those configured for ngs_mapping
        tools_somatic_variant_calling: []  # default to those configured for somatic_variant_calling
        target_regions: # REQUIRED
        padding: 0  # Used for count the number of variants outside of exom + padding
        AF_ID: 'AF' # REQUIRED ID of allele frequency field used in vcf file
        ignore_regions: "" # hard mapping regions
        minimal_support_read: 1 # threshold for defining a variant that has minimal support reads
        limited_support_read: 5 # threshold for defining a variant that has limited support reads
"""


class SomaticVariantQCStepPart(BaseStepPart):
    """"""

    name = "SVQC_step"

    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)

    @dictify
    def get_input_files(self, action):
        self._validate_action(action)
        tpl = (
            "output/{mapper}.{var_caller}.{tumor_library}/out/"
            "{mapper}.{var_caller}.{tumor_library}"
        )
        key_ext = {
            "full_vcf": ".full.vcf.gz",
            "full_vcf_tbi": ".full.vcf.gz.tbi",
            "passed_vcf": ".vcf.gz",
            "passed_vcf_tbi": ".vcf.gz.tbi",
        }
        variant_calling = self.parent.sub_workflows["somatic_variant_calling"]
        for key, ext in key_ext.items():
            yield key, variant_calling(tpl + ext)

    @dictify
    def get_output_files(self, action):
        # Validate action
        self._validate_action(action)
        prefix = (
            "work/{mapper}.{var_caller}.variantsqc.{tumor_library}/out/"
            "{mapper}.{var_caller}.variantsqc.{tumor_library}"
        )
        key_ext = {"json": ".json"}
        for key, ext in key_ext.items():
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action=action)
        prefix = (
            "work/{mapper}.{var_caller}.variantsqc.{tumor_library}/log/"
            "{mapper}.{var_caller}.variantsqc.{tumor_library}"
        )

        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_resource_usage(self, action):
        """Get Resource Usage
        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str
        :return: Returns ResourceUsage for step.
        :raises UnsupportedActionException: if action not in class defined list of valid actions.
        """
        self._validate_action(action=action)
        mem_mb = 4 * 1024  # 4GB
        return ResourceUsage(
            threads=2,
            time="1:00:00",  # 1 hour
            memory=f"{mem_mb}M",
        )


class SomaticVariantQCWorkflow(BaseStep):
    """Perform gathering variant information"""

    name = "somaticvariantchecking"
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
            (SomaticVariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((SomaticVariantQCStepPart, LinkOutStepPart))
        # Register sub workflows
        self.register_sub_workflow(
            "somatic_variant_calling",
            self.w_config["step_config"]["somatic_variant_checking"][
                "path_somatic_variant_calling"
            ],
        )
        # Copy over "tools" setting from somatic_variant_calling/ngs_mapping if not set here
        if not self.w_config["step_config"]["somatic_variant_checking"]["tools_ngs_mapping"]:
            self.w_config["step_config"]["somatic_variant_checking"][
                "tools_ngs_mapping"
            ] = self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"]
        if not self.w_config["step_config"]["somatic_variant_checking"][
            "tools_somatic_variant_calling"
        ]:
            self.w_config["step_config"]["somatic_variant_checking"][
                "tools_somatic_variant_calling"
            ] = self.w_config["step_config"]["somatic_variant_calling"]["tools"]

    @listify
    def get_result_files(self):
        callers = set(
            self.w_config["step_config"]["somatic_variant_checking"][
                "tools_somatic_variant_calling"
            ]
        )
        name_pattern = "{mapper}.{caller}.variantsqc.{tumor_library.name}"
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "out", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["somatic_variant_checking"]["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
            ext=EXT_VALUES,
        )
        yield from self._yield_result_files_matched(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            mapper=self.w_config["step_config"]["somatic_variant_checking"]["tools_ngs_mapping"],
            caller=callers & set(SOMATIC_VARIANT_CALLERS_MATCHED),
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
                if not sample_pair.tumor_sample.dna_ngs_library:
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

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "somatic_variant_checking", "path_somatic_variant_calling"),
            "Path to variant calling not configured but required for somatic variant qc",
        )

        self.ensure_w_config(
            ("step_config", "somatic_variant_checking", "target_regions"),
            "Path to target_regions file (bed format)"
            "not configured but required for somatic variant qc",
        )

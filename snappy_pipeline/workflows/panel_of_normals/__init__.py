# -*- coding: utf-8 -*-
"""Implementation of the ``panel_of_normals`` step

The ``panel_of_normals`` step takes as the input the results of the ``ngs_mapping`` step
(aligned reads in BAM format) and creates background information for somatic variant calling.
and/or somatic copy number calling. This background information is summarized as a
panel of normals.

Usually, the ``panel_of_normals`` step is required by somatic variant calling or
somatic copy number calling tools.

==========
Step Input
==========

The somatic variant calling step uses Snakemake sub workflows for using the result of the
``ngs_mapping`` step.

===========
Step Output
===========

For each panel of normals tool, the step outputs one set of files describing the panel.
For example, the ``mutect2`` panel of normal generates ``{mapper}.mutect2.pon.vcf.gz``
and associated files (md5 sums indices).

The normals that have been used, as well as the individual files (for example
vcf files for each normal) are kept in the ``work`` directory. This enables the
augmentation of the panel by new files when they become available.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_panel_of_normals.rst

=====================================
Panel of normals generation for tools
=====================================

- Panel of normal for ``mutect2`` somatic variant caller

=======
Reports
=======

Currently, no reports are generated.
"""

import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "tbi", "vcf_md5", "tbi_md5")

#: Names of the tools that might use panel of normals
TOOLS = ["mutect2"]

#: Default configuration for the somatic_variant_calling schema
DEFAULT_CONFIG = r"""
# Default configuration somatic_variant_calling
step_config:
  panel_of_normals:
    tools: ['mutect2']  # REQUIRED - available: 'mutect2'
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    ignore_chroms:             # patterns of chromosome names to ignore
    - NC_007605    # herpes virus
    - hs37d5       # GRCh37 decoy
    - chrEBV       # Eppstein-Barr Virus
    - '*_decoy'    # decoy contig
    - 'HLA-*'      # HLA genes
    - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
    # Configuration for mutect2
    mutect2:
      path_normals_list: null    # Optional file listing libraries to include in panel
      germline_resource: REQUIRED
      # Java options
      java_options: ' -Xmx16g '
      # Parallelization configuration
      num_cores: 2               # number of cores to use locally
      window_length: 100000000   # split input into windows of this size, each triggers a job
      num_jobs: 500              # number of windows to process in parallel
      use_profile: true          # use Snakemake profile for parallel processing
      restart_times: 5           # number of times to re-launch jobs in case of failure
      max_jobs_per_second: 2     # throttling of job creation
      max_status_checks_per_second: 10 # throttling of status checks
      debug_trunc_tokens: 0      # truncation to first N tokens (0 for none)
      keep_tmpdir: never         # keep temporary directory, {always, never, onerror}
      job_mult_memory: 1         # memory multiplier
      job_mult_time: 1           # running time multiplier
      merge_mult_memory: 1       # memory multiplier for merging
      merge_mult_time: 1         # running time multiplier for merging
"""


class PanelOfNormalsStepPart(BaseStepPart):
    """Base class for panel of normals step parts

    Two steps: the preparation is done on each normal samples separately, and the panel creation
    merges all the individual results in the the panel.
    """

    #: Class available actions
    actions = ("prepare_panel", "create_panel")

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.normal_libraries = self._get_normal_libraries()

    def _get_normal_libraries(self):
        for sheet in self.parent.shortcut_sheets:
            for donor in sheet.donors:
                for _, bio_sample in donor.bio_samples.items():
                    if bio_sample.is_tumor:
                        continue
                    for _, test_sample in bio_sample.test_samples.items():
                        extraction_type = test_sample.extra_infos.get("extractionType", "DNA")
                        if extraction_type.lower() == "dna":
                            for _, ngs_library in test_sample.ngs_libraries.items():
                                yield ngs_library.name


class Mutect2StepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "mutect2"

    #: Class resource usage dictionary. Key: action type (string); Value: resource (ResourceUsage).
    resource_usage_dict = {
        "prepare_panel": ResourceUsage(
            threads=2,
            time="3-00:00:00",  # 3 days
            memory="8G",
        ),
        "create_panel": ResourceUsage(
            threads=2,
            time="08:00:00",  # 8 hours
            memory="30G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)
        if self.config["mutect2"]["path_normals_list"]:
            self.normal_libraries = []
            with open(self.config["mutect2"]["path_normals_list"], "rt") as f:
                for line in f:
                    self.normal_libraries.append(line.strip())

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return self.resource_usage_dict.get(action)

    def get_input_files(self, action):
        """Return input files for mutect2 variant calling"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "prepare_panel": self._get_input_files_prepare_panel,
            "create_panel": self._get_input_files_create_panel,
        }
        return mapping[action]

    def _get_input_files_prepare_panel(self, wildcards):
        """Helper wrapper function for single sample panel preparation"""
        # Get shorcut to Snakemake sub workflow
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {"normal_bam": bam, "normal_bai": bam + ".bai"}

    def _get_input_files_create_panel(self, wildcards):
        """Helper wrapper function for merging individual results & panel creation"""
        paths = []
        tpl = "work/{mapper}.{tool}.prepare_panel/out/{normal_library}.vcf.gz"
        for normal in self.normal_libraries:
            paths.append(tpl.format(normal_library=normal, tool=self.name, **wildcards))
        return {"normals": paths}

    @dictify
    def get_output_files(self, action):
        """Return panel of normal files"""
        self._validate_action(action)
        ext_dict = {
            "vcf": "vcf.gz",
            "vcf_md5": "vcf.gz.md5",
            "tbi": "vcf.gz.tbi",
            "tbi_md5": "vcf.gz.tbi.md5",
        }
        tpls = {
            "prepare_panel": "work/{{mapper}}.{tool}.prepare_panel/out/{{normal_library}}.{ext}",
            "create_panel": "work/{{mapper}}.{tool}.create_panel/out/{{mapper}}.{tool}.panel_of_normals.{ext}",
        }
        for key, ext in ext_dict.items():
            yield key, tpls[action].format(tool=self.name, ext=ext)

    @dictify
    def get_log_file(self, action):
        """Return panel of normal files"""
        self._validate_action(action)
        ext_dict = {
            "conda_list": "conda_list.txt",
            "conda_list_md5": "conda_list.txt.md5",
            "conda_info": "conda_info.txt",
            "conda_info_md5": "conda_info.txt.md5",
            "log": "log",
            "log_md5": "log.md5",
        }
        tpls = {
            "prepare_panel": "work/{{mapper}}.{tool}.prepare_panel/log/{{normal_library}}.{ext}",
            "create_panel": "work/{{mapper}}.{tool}.create_panel/log/{{mapper}}.{tool}.panel_of_normals.{ext}",
        }
        for key, ext in ext_dict.items():
            yield key, tpls[action].format(tool=self.name, ext=ext)


class PanelOfNormalsWorkflow(BaseStep):
    """Creates a panel of normals"""

    # Workflow name
    name = "panel_of_normals"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

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
            (NgsMappingWorkflow,),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                Mutect2StepPart,
                LinkOutStepPart,
            )
        )
        # Initialize sub-workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all bio samples in all sample sheets.
        """
        # Panel of normals
        files = []

        # yield from expand(
        files.extend(
            expand(
                os.path.join(
                    "output",
                    "{mapper}.{caller}.create_panel",
                    "out",
                    "{mapper}.{caller}.panel_of_normals" + "{ext}",
                ),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=set(self.config["tools"]) & set(TOOLS),
                ext=EXT_VALUES,
            )
        )
        # yield from expand(
        files.extend(
            expand(
                os.path.join(
                    "output",
                    "{mapper}.{caller}.create_panel",
                    "log",
                    "{mapper}.{caller}.panel_of_normals" + "{ext}",
                ),
                mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                caller=set(self.config["tools"]) & set(TOOLS),
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                ),
            )
        )
        return files

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "panel_of_normals", "path_ngs_mapping"),
            "Path to NGS mapping not configured but required for somatic variant calling",
        )

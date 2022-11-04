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

from collections import OrderedDict
import os
import random
import sys

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
    tools: ['mutect2']
    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    size: 10
    shuffle_seed: 1234567
    # Configuration for mutect2
    mutect2:
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
      ignore_chroms:             # patterns of chromosome names to ignore
      - NC_007605    # herpes virus
      - hs37d5       # GRCh37 decoy
      - chrEBV       # Eppstein-Barr Virus
      - '*_decoy'    # decoy contig
      - 'HLA-*'      # HLA genes
      - 'GL000220.*' # Contig with problematic, repetitive DNA in GRCh37
"""


class PanelOfNormalsStepPart(BaseStepPart):
    """Base class for panel of normals step parts

    Variant calling is performed on matched cancer bio sample pairs.  That is, the primary NGS
    library for the primary bio sample is used for each cancer bio sample (paired with the primary
    normal bio sample's primary NGS library).
    """

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from cancer bio sample name to matched cancer sample
        self.tumor_ngs_library_to_sample_pair = OrderedDict()
        for sheet in self.parent.shortcut_sheets:
            self.tumor_ngs_library_to_sample_pair.update(
                sheet.all_sample_pairs_by_tumor_dna_ngs_library
            )

    def get_input_files(self, action):
        def input_function(wildcards):
            """Helper wrapper function"""
            print("DEBUG- input function [Part] wildcards = {}".format(wildcards), file=sys.stderr)
            # Get shorcut to Snakemake sub workflow
            ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
            # Get names of primary libraries of the selected cancer bio sample and the
            # corresponding primary normal sample
            normal_base_path = (
                "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}".format(**wildcards)
            )
            return {
                "normal_bam": ngs_mapping(normal_base_path + ".bam"),
                "normal_bai": ngs_mapping(normal_base_path + ".bam.bai"),
            }

        assert action == "run", "Unsupported actions"
        return input_function

    def get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        pair = self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]
        return pair.normal_sample.dna_ngs_library.name

    def get_output_files(self, action):
        """Return output files that all somatic variant calling sub steps must
        return (VCF + TBI file)
        """
        assert action == "run"
        return dict(
            zip(EXT_NAMES, expand(self.base_path_out, var_caller=[self.name], ext=EXT_VALUES))
        )

    def get_log_file(self, action):
        raise NotImplementedError("Panel of normals log file generation not implemented")


def get_panel_of_normals(filename, sheets, size, seed):
    """Reads from the normals list file, create it if necessary"""

    # Create the file if missing
    if not os.path.exists(filename):
        libraries = []
        for sheet in sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    if not bio_sample.extra_infos["isTumor"]:
                        libraries.append(bio_sample.dna_ngs_library.name)
        libraries.sort()
        random.seed(seed)
        random.shuffle(libraries)
        print("DEBUG- get_panel_of_normals bam files = {}".format(libraries), file=sys.stderr)
        if not os.path.exists(os.path.dirname(filename)):
            os.mkdir(os.path.dirname(filename))
        f = open(filename, "x")
        for normal_library in libraries[:size]:
            print(normal_library, file=f)
        f.close()

    normals = []
    with open(filename, "r") as f:
        for line in f:
            normals.append(line.rstrip())
    return normals


class Mutect2StepPart(PanelOfNormalsStepPart):
    """Somatic variant calling with MuTect 2"""

    #: Step name
    name = "mutect2"

    #: Class available actions
    actions = ("run", "prepare_panel", "create_panel")

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
        # TODO: Value set to default, maybe insufficient.
        "run": ResourceUsage(
            threads=1,
            time="01:00:00",
            memory="2G",
        ),
    }

    def __init__(self, parent):
        super().__init__(parent)

    def check_config(self):
        if self.name not in self.config["tools"]:
            return  # Mutect not enabled, skip
        self.parent.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for %s" % (self.name,),
        )

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
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        tpl = "output/{mapper}.{normal_library}/out/{mapper}.{normal_library}.bam"
        bam = ngs_mapping(tpl.format(**wildcards))
        return {"normal_bam": [bam], "normal_bai": [bam + ".bai"]}

    def _get_input_files_create_panel(self, wildcards):
        normals = self._get_panel_of_normals(wildcards)
        vcf_paths = []
        tpl = "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz"
        for normal in normals:
            vcf_paths.append(tpl.format(normal_library=normal, **wildcards))
        return {"txt": self._get_output_files_select_panel()["normal_list"], "vcf": vcf_paths}

    def get_output_files(self, action):
        """Return output files for mutect2 panel of normal creation"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "prepare_panel": self._get_output_files_prepare_panel,
            "create_panel": self._get_output_files_create_panel,
        }
        return mapping[action]()

    @dictify
    def _get_output_files_select_panel(self):
        # TODO: remove methods associated with action 'select_panel' if not used.
        yield "normal_list", "work/{mapper}.mutect2.select_panel.txt"

    @staticmethod
    def _get_output_files_prepare_panel():
        # TODO: Potential extension error in output files, `vcf.tbi.gz` instead of `vcf.gz.tbi`.
        return {
            "vcf": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz",
            "vcf_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.md5",
            "tbi": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.tbi.gz",
            "tbi_md5": "work/{mapper}.mutect2.prepare_panel/out/{normal_library}.vcf.gz.tbi.md5",
        }

    @staticmethod
    def _get_output_files_create_panel():
        base_name = "work/{mapper}.mutect2.create_panel/out/{mapper}.mutect2.panel_of_normals"
        return {
            "vcf": base_name + ".vcf.gz",
            "vcf_md5": base_name + ".vcf.gz.md5",
            "tbi": base_name + ".vcf.gz.tbi",
            "tbi_md5": base_name + ".vcf.gz.tbi.md5",
        }

    def get_log_file(self, action):
        """Return log files for mutect2 panel of normal creation"""
        # Validate action
        self._validate_action(action)
        mapping = {
            "prepare_panel": self._get_log_file_prepare_panel,
            "create_panel": self._get_log_file_create_panel,
        }
        return mapping[action]()

    @dictify
    def _get_log_file_prepare_panel(self):
        base_name = "work/{mapper}.mutect2.prepare_panel/log/{mapper}.mutect2.{normal_library}"
        return {
            "log": base_name + ".log",
            "log_md5": base_name + ".log.md5",
            "conda_info": base_name + ".conda_info.txt",
            "conda_info_md5": base_name + ".conda_info.txt.md5",
            "conda_list": base_name + ".conda_list.txt",
            "conda_list_md5": base_name + ".conda_list.txt.md5",
        }

    @dictify
    def _get_log_file_create_panel(self):
        base_name = "work/{mapper}.mutect2.create_panel/log/{mapper}.mutect2.create_panel"
        return {
            "log": base_name + ".log",
            "log_md5": base_name + ".log.md5",
            "conda_info": base_name + ".conda_info.txt",
            "conda_info_md5": base_name + ".conda_info.txt.md5",
            "conda_list": base_name + ".conda_list.txt",
            "conda_list_md5": base_name + ".conda_list.txt.md5",
        }

    @listify
    def _get_panel_of_normals(self, wildcards):
        """Return list of bam files in the panel of normals"""
        return get_panel_of_normals(
            "work/{mapper}.mutect2.select_panel.txt".format(**wildcards),
            self.parent.shortcut_sheets,
            self.config["size"],
            self.config["shuffle_seed"],
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
                    "{mapper}.{caller}.create_panel" + "{ext}",
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

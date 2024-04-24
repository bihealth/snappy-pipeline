# -*- coding: utf-8 -*-
"""Implementation of the ``variant_annotation`` step

The ``variant_annotation`` step takes as the input the results of the ``variant_calling`` step
(called germline variants in vcf.gz format) and annotates the variants, e.g., using VEP.

==========
Stability
==========

TBD

==========
Step Input
==========

The variant annotation step uses Snakemake sub workflows for using the result of the
``variant_calling`` step.

===========
Step Output
===========

TBD

====================
Global Configuration
====================

TBD

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_variant_annotation.rst

============================
Available Variant Annotators
============================

The following variant annotator is currently available:

- ``"vep"`` : See the `software documentation
  <https://www.ensembl.org/info/docs/tools/vep/script/index.html>`_ for more details

=======
Reports
=======

N/A
"""

import re
from itertools import chain

from biomedsheets.shortcuts import GermlineCaseSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.common import SnakemakeListItemsGenerator
from snappy_pipeline.workflows.abstract.exceptions import InvalidConfigurationException
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.variant_calling import GetResultFilesMixin, VariantCallingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Valid tools for variant annotation.
VARIANT_ANNOTATORS = ("vep",)

#: Extensions of files to create as main payload
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

# TODO: the number of restart times is high because tabix in HTSJDK/Jannovar is flaky...

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = r"""
# Default configuration variant_annotation
step_config:
  variant_annotation:
    path_variant_calling: ../variant_calling
    tools:
      - vep
    vep:
      # We will always run VEP in cache mode.  You have to provide the directory to the
      # cache to use (VEP would be ``~/.vep``).
      cache_dir: null # OPTIONAL
      # The cache version to use.  gnomAD v2 used 85, gnomAD v3.1 uses 101.
      cache_version: "85"
      # The assembly to use.  gnomAD v2 used "GRCh37", gnomAD v3.1 uses "GRCh38".
      assembly: "GRCh37"
      # The flag selecting the transcripts.  One of "gencode_basic", "refseq", and "merged".
      tx_flag: "gencode_basic"
      # Number of threads to use with forking, set to 0 to disable forking.
      num_threads: 16
      # Additional flags.
      more_flags: "--af_gnomade --af_gnomadg"
      # The --buffer_size parameter
      buffer_size: 100000
"""


class VepStepPart(GetResultFilesMixin, BaseStepPart):
    """Annotate VCF files using VEP"""

    #: Step name
    name = "vep"

    #: Class available actions
    actions = ("run",)

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        self._validate_action(action)
        token = "{mapper}.{var_caller}.{library_name}"
        variant_calling = self.parent.modules["variant_calling"]
        return {
            "vcf": variant_calling(f"output/{token}/out/{token}.vcf.gz"),
            "vcf_tbi": variant_calling(f"output/{token}/out/{token}.vcf.gz.tbi"),
        }

    @dictify
    def get_output_files(self, action):
        """Return output files for the filtration"""
        self._validate_action(action)
        token = "{mapper}.{var_caller}.vep.{library_name}"
        work_files = {
            "vcf": f"work/{token}/out/{token}.vcf.gz",
            "vcf_md5": f"work/{token}/out/{token}.vcf.gz.md5",
            "vcf_tbi": f"work/{token}/out/{token}.vcf.gz.tbi",
            "vcf_tbi_md5": f"work/{token}/out/{token}.vcf.gz.tbi.md5",
        }
        yield from work_files.items()
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(work_files.values(), self.get_log_file("run").values())
            ],
        )

    def get_extra_kv_pairs(self):
        return {"var_caller": self.parent.w_config["step_config"]["variant_calling"]["tools"]}

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        token = "{mapper}.{var_caller}.vep.{library_name}"
        prefix = f"work/{token}/log/{token}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, f"{prefix}{ext}"
            yield f"{key}_md5", f"{prefix}{ext}.md5"

    def get_resource_usage(self, action) -> ResourceUsage:
        self._validate_action(action)
        num_threads = self.config[self.name]["num_threads"]
        return ResourceUsage(
            threads=num_threads,
            time="1-00",
            memory=f"{2 * num_threads}G",
        )


class VariantAnnotationWorkflow(BaseStep):
    """Perform germline variant annotation"""

    name = "variant_annotation"
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
            (VariantCallingWorkflow, NgsMappingWorkflow),
        )
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((VepStepPart,))
        # Register sub workflows
        self.register_module(
            "ngs_mapping", self.w_config["step_config"]["variant_calling"]["path_ngs_mapping"]
        )
        self.register_module("variant_calling", self.config["path_variant_calling"])

    @listify
    def get_result_files(self) -> SnakemakeListItemsGenerator:
        for tool in self.config["tools"]:
            yield from self.sub_steps[tool].get_result_files()

    def check_config(self):
        """Check that the path to the NGS mapping is present"""
        self.ensure_w_config(
            ("step_config", "variant_annotation", "path_variant_calling"),
            "Path to variant calling not configured but required for variant annotation",
        )
        self.ensure_w_config(
            ("static_data_config", "reference", "path"),
            "Path to reference FASTA not configured but required for variant calling",
        )
        # Check that only valid tools are selected
        selected = set(self.w_config["step_config"]["variant_annotation"]["tools"])
        invalid = list(sorted(selected - set(VARIANT_ANNOTATORS)))
        if invalid:
            raise InvalidConfigurationException(f"Invalid variant callers selected: {invalid}")

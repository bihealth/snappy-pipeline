# -*- coding: utf-8 -*-
"""Implementation of the ``somatic_gene_fusion_calling`` step

The somatic_gene_fusion calling step allows for the detection of gene fusions from RNA-seq data
in cancer.  The wrapped tools start at the raw RNA-seq reads and generate filtered lists of
predicted gene fusions.

==========
Step Input
==========

Gene fusion calling starts at the raw RNA-seq reads.  Thus, the input is very similar to one of
:ref:`ngs_mapping step <step_ngs_mapping>`.

See :ref:`ngs_mapping_step_input` for more information.

===========
Step Output
===========

.. note:: TODO

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_gene_fusion_calling.rst

=============================
Available Gene Fusion Callers
=============================

- ``fusioncatcher``

"""

import os

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInPathGenerator,
    LinkInStep,
    LinkOutStepPart,
    ResourceUsage,
    get_ngs_library_folder_name,
)

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>"

#: HLA typing tools
GENE_FUSION_CALLERS = ("fusioncatcher", "jaffa", "pizzly", "hera", "star_fusion")

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = r"""
step_config:
  somatic_gene_fusion_calling:
    tools: ['fusioncatcher', 'jaffa']
    fusioncatcher:
      data_dir: REQUIRED   # REQUIRED
      configuration: null  # optional
      num_threads: 16
    pizzly:
      kallisto_index: REQUIRED    # REQUIRED
      transcripts_fasta: REQUIRED # REQUIRED
      annotations_gtf: REQUIRED       # REQUIRED
      kmer_size: 31
    hera:
      path_index: REQUIRED   # REQUIRED
      path_genome: REQUIRED  # REQUIRED
    star_fusion:
      path_ctat_resource_lib: REQUIRED
    defuse:
      path_dataset_directory: REQUIRED
""".lstrip()


class SomaticGeneFusionCallingStepPart(BaseStepPart):
    """Base class for somatic gene fusion calling"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_in = "work/input_links/{library_name}"
        self.base_path_out = "work/{name}.{{library_name}}/out/.done".format(name=self.name)
        # Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_set_infos, self.parent.config_lookup_paths
        )

    @dictify
    def get_input_files(self, action):
        """Return input files"""
        # Validate action
        self._validate_action(action)
        yield "done", os.path.join(self.base_path_in, ".done")

    @dictify
    def get_output_files(self, action):
        """Return output files that all read mapping sub steps must return (BAM + BAI file)"""
        # Validate action
        self._validate_action(action)
        yield "done", touch(self.base_path_out)

    def get_log_file(self, action):
        """Return path to log file"""
        # Validate action
        self._validate_action(action)
        return "work/{name}.{{library_name}}/log/snakemake.gene_fusion_calling.log".format(
            name=self.name
        )

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        _ = library_name
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        pattern_set_keys = ("right",) if prefix.startswith("right-") else ("left",)
        for _, path_infix, filename in self.path_gen.run(folder_name, pattern_set_keys):
            yield os.path.join(self.base_path_in, path_infix, filename).format(**wildcards)


class FusioncatcherStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Fusioncatcher"""

    #: Step name
    name = "fusioncatcher"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            tumor = flatten(
                zip(
                    sorted(self._collect_reads(wildcards, wildcards.library_name, "")),
                    sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")),
                )
            )
            tumor = list(map(os.path.abspath, tumor))
            return {"normal": [], "tumor": tumor}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{7500 * 4}M",
        )


class JaffaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using JAFFA"""

    #: Step name
    name = "jaffa"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{40 * 1024 * 4}M",
        )


class PizzlyStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Kallisto+Pizzly"""

    #: Step name
    name = "pizzly"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{20 * 1024 * 4}M",
        )


class StarFusionStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using STAR-Fusion"""

    #: Step name
    name = "star_fusion"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=4,
            time="5-00:00:00",  # 5 days
            memory=f"{30 * 1024 * 4}M",
        )


class DefuseStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Defuse"""

    #: Step name
    name = "defuse"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="5-00:00:00",  # 5 days
            memory=f"{10 * 1024 * 8}M",
        )


class HeraStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Hera"""

    #: Step name
    name = "hera"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"left": left, "right": right}

        assert action == "run", "Unsupported actions"
        return args_function

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=8,
            time="5-00:00:00",  # 5 days
            memory=f"{20 * 1024 * 8}M",
        )


class SomaticGeneFusionCallingWorkflow(BaseStep):
    """Perform somatic gene fusion calling"""

    #: Workflow name
    name = "somatic_gene_fusion_calling"

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific
        one
        """
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(workflow, config, config_lookup_paths, config_paths, workdir)
        self.register_sub_step_classes(
            (
                FusioncatcherStepPart,
                JaffaStepPart,
                PizzlyStepPart,
                HeraStepPart,
                StarFusionStepPart,
                DefuseStepPart,
                LinkInStep,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        name_pattern = "{fusion_caller}.{ngs_library.name}"
        for fusion_caller in self.config["tools"]:
            for sheet in self.shortcut_sheets:
                for donor in sheet.donors:
                    for bio_sample in donor.bio_samples.values():
                        for _test_sample in bio_sample.test_samples.values():
                            ngs_library = bio_sample.rna_ngs_library
                            if ngs_library is None:
                                break
                            name_pattern_value = name_pattern.format(
                                fusion_caller=fusion_caller, ngs_library=ngs_library
                            )
                            yield os.path.join("output", name_pattern_value, "out", ".done")

    def check_config(self):
        """Check that the required configurations are present."""
        # TODO: implement check for REQUIRED configurations.
        pass

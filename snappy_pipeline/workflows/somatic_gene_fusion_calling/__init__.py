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

.. note::

    The step requires a ``cancer_matched`` configuration & samplesheet files.
    This is an unnecessary requirement, which might be dropped in the future.

===========
Step Output
===========

There is no standard for reporting gene fusions, and therefore the output is different for all implemented tools.

``arriba`` returns two tab-separated files: ``arriba.<library name>.fusions.tsv`` & ``arriba.<library name>.discarded_fusions.tsv.gz``.
Both files list the affected genes, reads supporting the fusion & a confidence level.
Obviously, the discarded fusion file contains all hints of fusion that have been discarded because of insufficient evidence.

=====================
Default Configuration
=====================

The default configuration is as follows.

.. include:: DEFAULT_CONFIG_somatic_gene_fusion_calling.rst

=============================
Available Gene Fusion Callers
=============================

- ``arriba``
- ``fusioncatcher`` implementation is broken. The tool's computational resources requirements are so enormous that it might not be advisable to try re-enable it.
- the status of ``defuse``, ``hera``, ``jaffa`` & ``pizzly`` is unknown, they are probably currently broken or not implemented.
- the status of ``star_fusion`` is also unknown, but it apparently returns results fairly similar to ``arriba``, but not quite as accurate. ``arriba`` should be preferred.

"""

import os

from snakemake.io import touch

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
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

from .model import SomaticGeneFusionCalling as SomaticGeneFusionCallingConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: HLA typing tools
GENE_FUSION_CALLERS = (
    "arriba",
    "defuse",
    "fusioncatcher",
    "hera",
    "jaffa",
    "pizzly",
    "star_fusion",
)

#: Default configuration for the somatic_gene_fusion_calling step
DEFAULT_CONFIG = SomaticGeneFusionCallingConfigModel.default_config_yaml_string()


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
            self.parent.work_dir,
            self.parent.data_set_infos,
            self.parent.config_lookup_paths,
            preprocessed_path=self.config.path_link_in,
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
        folder_name = get_ngs_library_folder_name(self.parent.sheets, wildcards.library_name)
        if self.config.path_link_in:
            folder_name = library_name
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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


class ArribaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using arriba"""

    #: Step name
    name = "arriba"

    def get_args(self, action):
        """Return function that maps wildcards to dict for input files"""

        def flatten(lst):
            return [x for pair in lst for x in pair]

        def args_function(wildcards):
            # TODO: wildcards.library_name is tumor_library_name
            left = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "")))
            right = list(sorted(self._collect_reads(wildcards, wildcards.library_name, "right-")))
            return {"input": {"reads_left": left, "reads_right": right}}

        assert action == "run", "Unsupported actions"
        return args_function

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        base_path_out = "work/{name}.{{library_name}}/out/{name}.{{library_name}}.{ext}"
        key_ext = (
            ("fusions", "fusions.tsv"),
            ("discarded", "discarded_fusions.tsv.gz"),
        )
        for key, ext in key_ext:
            yield key, base_path_out.format(name=self.name, ext=ext)
            yield key + "_md5", base_path_out.format(name=self.name, ext=ext) + ".md5"
        yield "done", "work/arriba.{library_name}/out/.done"

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = "work/{name}.{{library_name}}/log/{name}.{{library_name}}".format(name=self.name)
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"
        prefix = "work/{name}.{{library_name}}/log/".format(name=self.name)
        key_ext = (
            ("out", "Log.out"),
            ("final", "Log.final.out"),
            ("std", "Log.std.out"),
            ("SJ", "SJ.out.tab"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        """Get resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        # Validate action
        self._validate_action(action)
        return ResourceUsage(
            threads=self.config.arriba.num_threads, time="24:00:00", memory=f"{96 * 1024}M"
        )  # 1 day


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
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=SomaticGeneFusionCallingConfigModel,
        )
        self.register_sub_step_classes(
            (
                FusioncatcherStepPart,
                JaffaStepPart,
                PizzlyStepPart,
                HeraStepPart,
                StarFusionStepPart,
                DefuseStepPart,
                ArribaStepPart,
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
        # Convert sheet parsing into method
        library_names_list = list(self._get_all_rna_ngs_libraries())
        # Get results
        name_pattern = "{fusion_caller}.{ngs_library}"
        for fusion_caller in self.config.tools:
            for ngs_library in library_names_list:
                # Constant to all callers
                name_pattern_value = name_pattern.format(
                    fusion_caller=fusion_caller, ngs_library=ngs_library
                )
                yield os.path.join("output", name_pattern_value, "out", ".done")
                # Caller specific stuff...
                if fusion_caller == "arriba":
                    yield from self._yield_arriba_files(ngs_library)
                else:
                    yield os.path.join(
                        "output", name_pattern_value, "log", "snakemake.gene_fusion_calling.log"
                    )

    def _get_all_rna_ngs_libraries(self):
        for sheet in self.shortcut_sheets:
            for donor in sheet.donors:
                for _, bio_sample in donor.bio_samples.items():
                    for _, test_sample in bio_sample.test_samples.items():
                        extraction_type = test_sample.extra_infos.get("extractionType", "DNA")
                        if extraction_type.lower() == "rna":
                            for _, ngs_library in test_sample.ngs_libraries.items():
                                yield ngs_library.name

    def _yield_arriba_files(self, ngs_library):
        tpl = "output/arriba.{library_name}/out/arriba.{library_name}.{ext}"
        for ext in ("fusions.tsv", "discarded_fusions.tsv.gz"):
            yield tpl.format(library_name=ngs_library, ext=ext)
            yield tpl.format(library_name=ngs_library, ext=ext + ".md5")
        tpl = "output/arriba.{library_name}/log/arriba.{library_name}.{ext}"
        for ext in ("log", "conda_list.txt", "conda_info.txt"):
            yield tpl.format(library_name=ngs_library, ext=ext)
            yield tpl.format(library_name=ngs_library, ext=ext + ".md5")
        tpl = "output/arriba.{library_name}/log/{ext}"
        for ext in ("Log.out", "Log.std.out", "Log.final.out", "SJ.out.tab"):
            yield tpl.format(library_name=ngs_library, ext=ext)
            yield tpl.format(library_name=ngs_library, ext=ext + ".md5")

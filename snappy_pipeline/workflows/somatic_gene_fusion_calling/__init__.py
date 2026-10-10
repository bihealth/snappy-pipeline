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

``arriba`` returns two tab-separated files: ``<library name>.fusions.tsv`` & ``<library name>.discarded_fusions.tsv.gz``.
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

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import SomaticGeneFusionCalling as SomaticGeneFusionCallingConfigModel
from .model import Tool

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


class SomaticGeneFusionCallingStepPart(BaseStepPart):
    """Base class for somatic gene fusion calling"""

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)
        self.base_path_out = "work/{library_name}/out/.done"

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return the ordered mates, and the upstream ``.done`` file when a task wrote the reads"""
        library_name = wildcards.library_name
        left = sorted(self._collect_reads(wildcards, library_name, ""))
        right = sorted(self._collect_reads(wildcards, library_name, "right-"))
        yield "reads_left", left
        if right:
            yield "reads_right", right
        written = {*left, *right}
        if done := [path for path in self.parent.reads_input(library_name) if path not in written]:
            yield "reads_done", done

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
        return "work/{library_name}/log/snakemake.gene_fusion_calling.log"

    def _collect_reads(self, wildcards, library_name, prefix):
        """Yield the path to reads

        Yields paths to right reads if prefix=='right-'
        """
        _ = wildcards
        mate = "right" if prefix.startswith("right-") else "left"
        for group in self.parent.read_groups(library_name):
            if mate in group.paths:
                yield group.paths[mate]


class FusioncatcherStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Fusioncatcher"""

    #: Step name
    name = "fusioncatcher"

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
            runtime="5d",  # 5 days
            mem=f"{7500 * 4}MB",
        )


class JaffaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using JAFFA"""

    #: Step name
    name = "jaffa"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "reference_files", self.config.jaffa.path_reference_files

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
            runtime="5d",  # 5 days
            mem=f"{40 * 1024 * 4}MB",
        )


class PizzlyStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Kallisto+Pizzly"""

    #: Step name
    name = "pizzly"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "kallisto_index", self.config.pizzly.kallisto_index
        yield "transcripts_fasta", self.config.pizzly.transcripts_fasta
        yield "annotations_gtf", self.config.pizzly.annotations_gtf

    def get_params(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            return {"kmer_size": self.config.pizzly.kmer_size}

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
            runtime="5d",  # 5 days
            mem=f"{20 * 1024 * 4}MB",
        )


class StarFusionStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using STAR-Fusion"""

    #: Step name
    name = "star_fusion"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "ctat_resource_lib", self.config.star_fusion.path_ctat_resource_lib

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
            runtime="5d",  # 5 days
            mem=f"{30 * 1024 * 4}MB",
        )


class DefuseStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Defuse"""

    #: Step name
    name = "defuse"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "dataset_directory", self.config.defuse.path_dataset_directory

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
            runtime="5d",  # 5 days
            mem=f"{10 * 1024 * 8}MB",
        )


class HeraStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Hera"""

    #: Step name
    name = "hera"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "genome", self.config.hera.path_genome
        yield "index", self.config.hera.path_index

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
            runtime="5d",  # 5 days
            mem=f"{20 * 1024 * 8}MB",
        )


class ArribaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using arriba"""

    #: Step name
    name = "arriba"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "reference", self.parent.get_upstream_paths("reference").fasta
        yield "features", self.parent.get_upstream_paths("features").gtf
        yield "index", self.config.arriba.path_index
        for key in ("blacklist", "known_fusions", "tags", "structural_variants", "protein_domains"):
            if path := getattr(self.config.arriba, key):
                yield key, path

    def get_params(self, action):
        """Return function that maps wildcards to dict for input files"""

        def args_function(wildcards):
            return {
                "trim_adapters": self.config.arriba.trim_adapters,
                "num_threads_trimming": self.config.arriba.num_threads_trimming,
                "num_threads": self.config.arriba.num_threads,
                "star_parameters": self.config.arriba.star_parameters,
            }

        assert action == "run", "Unsupported actions"
        return args_function

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        base_path_out = "work/{{library_name}}/out/{{library_name}}.{ext}"
        key_ext = (
            ("fusions", "fusions.tsv"),
            ("discarded", "discarded_fusions.tsv.gz"),
        )
        for key, ext in key_ext:
            path = base_path_out.format(ext=ext)
            yield key, path
            yield key + "_md5", path + ".md5"
        yield "done", "work/{library_name}/out/.done"

    @dictify
    def get_log_file(self, action):
        """Return dict of log files."""
        _ = action
        prefix = "work/{library_name}/log/{library_name}"
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"
        prefix = "work/{library_name}/log/"
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
            threads=self.config.arriba.num_threads, runtime="24h", mem=f"{96 * 1024}MB"
        )  # 1 day


class SomaticGeneFusionCallingWorkflow(BaseStep):
    """Perform somatic gene fusion calling"""

    #: Workflow name
    name = "somatic_gene_fusion_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "fusion", "rna"}))]

    config_model_class = SomaticGeneFusionCallingConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local fusion-calling output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"output/{lib}/out/.done"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        selected_tool = self.config.tool
        match selected_tool:
            case Tool.fusioncatcher:
                selected_sub_step = FusioncatcherStepPart
            case Tool.jaffa:
                selected_sub_step = JaffaStepPart
            case Tool.pizzly:
                selected_sub_step = PizzlyStepPart
            case Tool.hera:
                selected_sub_step = HeraStepPart
            case Tool.star_fusion:
                selected_sub_step = StarFusionStepPart
            case Tool.defuse:
                selected_sub_step = DefuseStepPart
            case Tool.arriba:
                selected_sub_step = ArribaStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {selected_tool}")
        self.register_sub_step_classes(
            (
                selected_sub_step,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        df = self.build_library_dataframe()
        if df.empty:
            return []

        fusion_tool = str(self.config.tool)
        name_pattern = "{library_name}"
        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type.lower() != "rna":
                continue
            name_pattern_value = name_pattern.format(library_name=library_name)
            yield os.path.join("output", name_pattern_value, "out", ".done")
            if fusion_tool == "arriba":
                yield from self._yield_arriba_files(library_name)
            else:
                yield os.path.join(
                    "output", name_pattern_value, "log", "snakemake.gene_fusion_calling.log"
                )

    def _yield_arriba_files(self, library_name):
        tpl = "output/{library_name}/out/{library_name}.{ext}"
        for ext in ("fusions.tsv", "discarded_fusions.tsv.gz"):
            yield tpl.format(library_name=library_name, ext=ext)
            yield tpl.format(library_name=library_name, ext=ext + ".md5")
        tpl = "output/{library_name}/log/{library_name}.{ext}"
        for ext in ("log", "conda_list.txt", "conda_info.txt"):
            yield tpl.format(library_name=library_name, ext=ext)
            yield tpl.format(library_name=library_name, ext=ext + ".md5")
        tpl = "output/{library_name}/log/{ext}"
        for ext in ("Log.out", "Log.std.out", "Log.final.out", "SJ.out.tab"):
            yield tpl.format(library_name=library_name, ext=ext)
            yield tpl.format(library_name=library_name, ext=ext + ".md5")

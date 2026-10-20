# -*- coding: utf-8 -*-
"""Somatic gene fusion calling from RNA-seq data

The step detects gene fusions in cancer RNA-seq data; a task takes all RNA libraries unless it sets ``library_selection``.
Each tool starts at the raw RNA-seq reads and predicts the fusions itself.
The reads come from the task named in ``depends_on.reads``.

The ``tool`` field selects the caller:

- ``arriba`` is the recommended caller.
  It needs the STAR index of ``depends_on.index`` (``{type: index, tags: [star, rna]}``, preferably STAR 2.7.10 or later).
  The optional ``depends_on.database`` (``{type: database, tags: [arriba]}``) may give the keys ``blacklist``, ``known_fusions``, ``tags``, ``structural_variants`` and ``protein_domains`` of the arriba release.
- ``star_fusion`` returns results similar to ``arriba``, but less accurate.
  It needs ``depends_on.database`` (``{type: database, tags: [star_fusion]}``) with the key ``ctat_resource_lib``.
- ``fusioncatcher`` is broken.
  Its resource requirements are so high that re-enabling it is probably not worth it.
- ``defuse``, ``hera``, ``jaffa`` and ``pizzly`` are not tested and probably broken.

The task gives one result per RNA library.
The sample sheet must be a cancer sheet.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/`` as symlinks to ``work/``.

- ``{library_name}/out/.done``: marks a finished library.
  It is the only file of every tool except ``arriba``.
- ``{library_name}/out/{library_name}.fusions.tsv`` and ``{library_name}/out/{library_name}.discarded_fusions.tsv.gz`` (``arriba``), each with a ``.md5`` file.
  The files list the affected genes, the supporting reads and a confidence level.
  The discarded fusions have insufficient evidence.

Downstream tasks read the key ``done``.
There is no standard format for gene fusions, so the other files depend on the tool.
"""

from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions
from snakemake.io import touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
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

    def _get_output_files_run(self):
        return {"done": touch("work/{library_name}/out/.done")}

    def _get_log_file_run(self):
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

    default_resource_usage = ResourceUsage(threads=4, runtime="5d", mem=f"{7500 * 4}MB")


class JaffaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using JAFFA"""

    #: Step name
    name = "jaffa"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "reference_files", self.config.jaffa.path_reference_files

    default_resource_usage = ResourceUsage(threads=4, runtime="5d", mem=f"{40 * 1024 * 4}MB")


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

    default_resource_usage = ResourceUsage(threads=4, runtime="5d", mem=f"{20 * 1024 * 4}MB")


class StarFusionStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using STAR-Fusion"""

    #: Step name
    name = "star_fusion"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "ctat_resource_lib", self.parent.get_upstream_paths("database").ctat_resource_lib

    default_resource_usage = ResourceUsage(threads=4, runtime="5d", mem=f"{30 * 1024 * 4}MB")


class DefuseStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Defuse"""

    #: Step name
    name = "defuse"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "dataset_directory", self.config.defuse.path_dataset_directory

    default_resource_usage = ResourceUsage(threads=8, runtime="5d", mem=f"{10 * 1024 * 8}MB")


class HeraStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using Hera"""

    #: Step name
    name = "hera"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "genome", self.config.hera.path_genome
        yield "index", self.config.hera.path_index

    default_resource_usage = ResourceUsage(threads=8, runtime="5d", mem=f"{20 * 1024 * 8}MB")


class ArribaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using arriba"""

    #: Step name
    name = "arriba"

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield from self.parent.reference_inputs().items()
        yield "features", self.parent.get_upstream_paths("features").gtf
        yield "index", self.parent.get_upstream_paths("index").index
        if self.config.depends_on.database:
            database = self.parent.get_upstream_paths("database")
            for key in type(database).model_fields:
                if path := getattr(database, key):
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
    def _get_output_files_run(self):
        key_ext = (
            ("fusions", "fusions.tsv"),
            ("discarded", "discarded_fusions.tsv.gz"),
        )
        for key, ext in key_ext:
            path = f"work/{{library_name}}/out/{{library_name}}.{ext}"
            yield key, path
            yield key + "_md5", path + ".md5"
        yield "done", "work/{library_name}/out/.done"

    @dictify
    def _get_log_file_run(self):
        """Return dict of log files."""
        yield "log", "work/{library_name}/log/{library_name}.log"
        prefix = "work/{library_name}/log/"
        key_ext = (
            ("out", "Log.out"),
            ("final", "Log.final.out"),
            ("std", "Log.std.out"),
            ("SJ", "SJ.out.tab"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(
            threads=self.config.arriba.num_threads, runtime="24h", mem=f"{96 * 1024}MB"
        )


class SomaticGeneFusionCallingWorkflow(BaseStep):
    """Perform somatic gene fusion calling"""

    #: Workflow name
    name = "somatic_gene_fusion_calling"
    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "fusion", "rna"}))]

    config_model_class = SomaticGeneFusionCallingConfigModel

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    #: RNA libraries, unless the task sets ``library_selection``
    library_default = "extraction_type == 'rna'"

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local fusion-calling output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"done": f"work/{lib}/out/.done"}

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
        self.register_sub_step_classes((selected_sub_step,))

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample
        sheets.
        """
        df = self.build_library_dataframe()
        if df.empty:
            return []

        for library_name in self.output_entities:
            row = df[df["library_name"] == library_name]
            if row.empty:
                continue
            extraction_type = row.iloc[0].get("extraction_type", "")
            if extraction_type.lower() != "rna":
                continue
            yield f"output/{library_name}/out/.done"
            if self.config.tool == Tool.arriba:
                for ext in ("fusions.tsv", "discarded_fusions.tsv.gz"):
                    yield f"output/{library_name}/out/{library_name}.{ext}"
                    yield f"output/{library_name}/out/{library_name}.{ext}.md5"

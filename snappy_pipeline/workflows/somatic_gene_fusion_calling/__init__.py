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

The task gives one result per RNA library.
The sample sheet must be a cancer sheet.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/`` as symlinks to ``work/``.

- ``arriba``: ``{library_name}/out/{library_name}.fusions.tsv`` and ``{library_name}/out/{library_name}.discarded_fusions.tsv.gz``, each with a ``.md5`` file.
  The files list the affected genes, the supporting reads and a confidence level.
  The discarded fusions have insufficient evidence.
  Downstream tasks read the key ``fusions``.
- ``star_fusion``: the directory ``{library_name}/out/``, which holds the STAR-Fusion output (``output/``) and its input reads (``input/``).
  Downstream tasks read the key ``out``.

There is no standard format for gene fusions, so the files depend on the tool.
"""

from snakemake.io import directory

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

from .model import SomaticGeneFusionCalling as SomaticGeneFusionCallingConfigModel
from .model import Tool

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: The output directory of a library
_OUT = "work/{library}/out"

#: The tables of arriba by output key: the fusions and the discarded candidates
ARRIBA_TABLES = {"fusions": "fusions.tsv", "discarded": "discarded_fusions.tsv.gz"}


class SomaticGeneFusionCallingStepPart(BaseStepPart):
    """Base class for somatic gene fusion calling"""

    #: Class available actions
    actions = ("run",)

    def _get_input_files_run(self, wildcards):
        """Return the ordered mates, or the directory of the task that wrote them"""
        return reads_input_files(self.parent, wildcards.library_name)

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/snakemake.gene_fusion_calling.log"}


class StarFusionStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using STAR-Fusion"""

    #: Step name
    name = "star_fusion"

    wrappers = {"run": "star_fusion"}
    result_keys = ("out",)
    default_resource_usage = ResourceUsage(threads=4, runtime="5d", mem=f"{30 * 1024 * 4}MB")

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        yield "ctat_resource_lib", self.parent.get_upstream_paths("database").ctat_resource_lib

    def _get_output_files_run(self):
        # STAR-Fusion writes many files that the step does not name
        return {"out": directory(_OUT.format(library="{library_name}"))}

    def _get_params_run(self, wildcards):
        """The reads written by another task, and the tool's ``extra_arguments``, which
        ``with_extra_arguments`` adds"""
        return reads_params(self.parent, wildcards.library_name)


class ArribaStepPart(SomaticGeneFusionCallingStepPart):
    """Somatic gene fusion calling from RNA-seq reads using arriba"""

    #: Step name
    name = "arriba"

    wrappers = {"run": "arriba/run"}
    result_keys = ("fusions", "fusions_md5", "discarded", "discarded_md5")

    @property
    def default_resource_usage(self) -> ResourceUsage:
        return ResourceUsage(
            threads=self.config.arriba.num_threads, runtime="24h", mem=f"{96 * 1024}MB"
        )

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

    @dictify
    def _get_output_files_run(self):
        out = _OUT.format(library="{library_name}")
        for key, suffix in ARRIBA_TABLES.items():
            path = f"{out}/{{library_name}}.{suffix}"
            yield key, path
            yield key + "_md5", path + ".md5"

    def _get_params_run(self, wildcards):
        return {
            "trim_adapters": self.config.arriba.trim_adapters,
            "num_threads_trimming": self.config.arriba.num_threads_trimming,
            "num_threads": self.config.arriba.num_threads,
            "star_parameters": self.config.arriba.star_parameters,
        } | reads_params(self.parent, wildcards.library_name)


class SomaticGeneFusionCallingWorkflow(BaseStep):
    """Perform somatic gene fusion calling"""

    #: Workflow name
    name = "somatic_gene_fusion_calling"

    produces = [DataSignature(DataType.VARIANTS, frozenset({"somatic", "fusion", "rna"}))]
    config_model_class = SomaticGeneFusionCallingConfigModel
    tool_parts = {Tool.star_fusion: StarFusionStepPart, Tool.arriba: ArribaStepPart}
    library_requirement = "kind == 'cancer' and extraction_type == 'rna'"
    library_default = ""

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local fusion-calling output paths for downstream consumers: the fusions table
        of arriba, the output directory of STAR-Fusion."""
        lib = kwargs.get("library_name", "{library_name}")
        out = _OUT.format(library=lib)
        if config.tool == Tool.arriba:
            return {"fusions": f"{out}/{lib}.{ARRIBA_TABLES['fusions']}"}
        return {"out": out}

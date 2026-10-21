# -*- coding: utf-8 -*-
"""Indices of a reference sequence

Builds the indices of a reference FASTA file for one mapper.
The ``tool`` field selects ``bwa``, ``bwa_mem2``, ``minimap2`` or ``star``.
STAR needs an RNA reference (``reference_molecule: rna``) and takes a gene annotation from ``depends_on.features`` to add splice junctions.
The other tools need a DNA reference.
The task has one result per project, not per library.

Outputs
-------

Below ``tasks/<task>/output/reference_index/out/``, with the prefix ``reference``:

- ``reference.fa.fai``, ``reference.dict`` and ``reference.fa.genome`` for every tool.
- ``bwa``: ``.amb``, ``.ann``, ``.bwt``, ``.pac`` and ``.sa``.
- ``bwa_mem2``: ``.0123``, ``.amb``, ``.ann``, ``.bwt.2bit.64`` and ``.pac``.
- ``minimap2``: ``.mmi``.
- ``star``: the directory ``reference.index/`` with the STAR genome files, plus the splice junction files if ``features`` is set.

Downstream tasks read the index through the key ``index``.
This is the prefix of the files for ``bwa`` and ``bwa_mem2``, the ``.mmi`` file for ``minimap2`` and the directory for ``star``.
The keys ``reference_fai``, ``reference_dict`` and ``reference_genome`` name the other three files.
"""

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ReferenceIndex as ReferenceIndexConfigModel
from .model import Tool

#: Files that ``STAR --runMode genomeGenerate`` writes into the index directory
STAR_INDEX_FILES = (
    "chrLength.txt",
    "chrName.txt",
    "chrNameLength.txt",
    "chrStart.txt",
    "Genome",
    "genomeParameters.txt",
    "SA",
    "SAindex",
)

#: Additional files with splice junctions from a gene annotation (``--sjdbGTFfile``)
STAR_SJDB_FILES = (
    "exonGeTrInfo.tab",
    "exonInfo.tab",
    "geneInfo.tab",
    "sjdbInfo.txt",
    "sjdbList.fromGTF.out.tab",
    "sjdbList.out.tab",
    "transcriptInfo.tab",
)


class BuildReferenceCommonStepPart(BaseStepPart):
    name = "common"
    actions = ("run",)

    default_resource_usage = ResourceUsage(threads=1, runtime="2h", mem="2GB")

    def _get_input_files_run(self, wildcards):
        return {"reference": self.parent.get_reference_path()}

    def _get_output_files_run(self):
        work_prefix = "work/reference_index/out/reference"
        return {
            "reference_fai": work_prefix + ".fa.fai",
            "reference_dict": work_prefix + ".dict",
            "reference_genome": work_prefix + ".fa.genome",
        }

    def _get_log_file_run(self):
        prefix = "work/reference_index/log/reference_index.common"
        return {"log": prefix + ".log"}


class _IndexToolStepPart(BaseStepPart):
    output_suffixes: tuple[str, ...] = ()

    def _get_input_files_run(self, wildcards):
        return {
            "reference": self.parent.get_reference_path(),
            "reference_fai": "work/reference_index/out/reference.fa.fai",
            "reference_dict": "work/reference_index/out/reference.dict",
            "reference_genome": "work/reference_index/out/reference.fa.genome",
        }

    def _get_output_files_run(self):
        work_prefix = "work/reference_index/out/reference"
        return {
            suffix.strip(".").replace(".", "_").replace("/", "_"): work_prefix + suffix
            for suffix in self.output_suffixes
        }

    def _get_log_file_run(self):
        prefix = "work/reference_index/log/reference_index.index"
        return {"log": prefix + ".log"}


class BwaIndexStepPart(_IndexToolStepPart):
    name = "bwa"
    actions = ("run",)
    output_suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")

    default_resource_usage = ResourceUsage(threads=8, runtime="24h", mem="16GB")

    @dictify
    def _get_params_run(self, wildcards):
        yield "algorithm", self.config.bwa.algorithm


class BwaMem2IndexStepPart(_IndexToolStepPart):
    name = "bwa_mem2"
    actions = ("run",)
    output_suffixes = (".0123", ".amb", ".ann", ".bwt.2bit.64", ".pac")

    default_resource_usage = ResourceUsage(threads=8, runtime="24h", mem="16GB")

    @dictify
    def _get_params_run(self, wildcards):
        yield "extra_args", " ".join(self.config.bwa_mem2.extra_args)


class Minimap2IndexStepPart(_IndexToolStepPart):
    name = "minimap2"
    actions = ("run",)
    output_suffixes = (".mmi",)

    default_resource_usage = ResourceUsage(threads=8, runtime="12h", mem="16GB")

    @dictify
    def _get_params_run(self, wildcards):
        yield "extra_args", " ".join(self.config.minimap2.extra_args)


class StarIndexStepPart(_IndexToolStepPart):
    name = "star"
    actions = ("run",)

    @property
    def output_suffixes(self) -> tuple[str, ...]:
        files = STAR_INDEX_FILES + (STAR_SJDB_FILES if self.config.depends_on.features else ())
        return tuple(f".index/{name}" for name in files)

    default_resource_usage = ResourceUsage(threads=16, runtime="24h", mem="64GB")

    @dictify
    def _get_params_run(self, wildcards):
        yield "extra_args", " ".join(self.config.star.extra_args)

    @dictify
    def _get_input_files_run(self, wildcards):
        yield from super()._get_input_files_run(wildcards).items()
        if self.config.depends_on.features:
            yield "features", self.parent.get_upstream_paths("features").gtf


class ReferenceIndexWorkflow(BaseStep):
    name = "reference_index"
    produces = [
        DataSignature(DataType.INDEX, frozenset({"bwa", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"bwa_mem2", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"minimap2", "dna"})),
        DataSignature(DataType.INDEX, frozenset({"star", "rna"})),
    ]

    config_model_class = ReferenceIndexConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        molecule = "rna" if config.tool == Tool.star else "dna"
        return (DataSignature(DataType.INDEX, frozenset({str(config.tool), molecule})),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        prefix = kwargs.get("prefix", "work/reference_index/out/reference")
        index = {Tool.minimap2: prefix + ".mmi", Tool.star: prefix + ".index"}
        return {
            "index": index.get(Tool(config.tool), prefix),
            "reference_fai": prefix + ".fa.fai",
            "reference_dict": prefix + ".dict",
            "reference_genome": prefix + ".fa.genome",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        tool_to_class = {
            Tool.bwa: BwaIndexStepPart,
            Tool.bwa_mem2: BwaMem2IndexStepPart,
            Tool.minimap2: Minimap2IndexStepPart,
            Tool.star: StarIndexStepPart,
        }
        selected_tool_class = tool_to_class[self.config.tool]
        self.register_sub_step_classes((BuildReferenceCommonStepPart, selected_tool_class))

    def get_reference_path(self) -> str:
        return self.get_upstream_paths("reference").fasta

    @listify
    def get_result_files(self):
        for tool in ("common", str(self.config.tool)):
            for path in self.get_output_files(tool, "run").values():
                yield path.replace("work/", "output/", 1)

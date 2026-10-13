# -*- coding: utf-8 -*-
"""Download of a reference sequence

Downloads one FASTA file from Ensembl, RefSeq or UCSC.
The ``source`` field selects the source.
The source section (``ensembl``, ``refseq`` or ``ucsc``) selects the species, build and file, or gives a ``url`` that replaces the generated URL.
``contigs`` and ``contigs_regex`` keep only some contigs of the file.
The task has one result per project, not per library.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``reference_download/out/reference.fa`` and ``reference_download/out/reference.fa.md5``

Downstream tasks read the FASTA file through the key ``fasta`` of their ``reference`` field.
The task produces data of type ``raw`` with the tags ``reference``, the source and the molecule (``dna``, ``rna`` or ``protein``).
"""

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ReferenceDownload as ReferenceDownloadConfigModel
from .model import Source


class _ReferenceDownloadStepPart(BaseStepPart):
    source: Source
    actions = ("run",)

    def _get_input_files_run(self, wildcards):
        return []

    default_resource_usage = ResourceUsage(threads=2, runtime="2h", mem="2GB")

    def _get_output_files_run(self):
        work_fasta = "work/reference_download/out/reference.fa"
        return {"fasta": work_fasta, "fasta_md5": work_fasta + ".md5"}

    def _base_args(self) -> dict:
        source_cfg = getattr(self.config, self.name)
        return {
            "source": self.name,
            "url": source_cfg.url,
            "species": source_cfg.species,
            "build": source_cfg.build,
            "assembly": source_cfg.assembly,
            "version": source_cfg.version,
            "contigs": source_cfg.contigs,
            "contigs_regex": source_cfg.contigs_regex,
            "molecule": source_cfg.molecule,
        }

    def _get_log_file_run(self):
        return {"log": f"work/reference_download/log/reference_download.{self.name}.log"}

    @dictify
    def _get_params_run(self, wildcards):
        yield from self._base_args().items()


class EnsemblReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ensembl"
    source = Source.ensembl

    @dictify
    def _get_params_run(self, wildcards):
        yield from self._base_args().items()
        cfg = self.config.ensembl
        yield "subset", cfg.subset
        yield "file_name", cfg.file_name
        yield "datatype", cfg.datatype
        yield "chromosome", cfg.chromosome
        yield "branch", cfg.branch


class RefseqReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "refseq"
    source = Source.refseq

    @dictify
    def _get_params_run(self, wildcards):
        yield from self._base_args().items()
        cfg = self.config.refseq
        yield "assembly_accession", cfg.assembly_accession
        yield "assembly_name", cfg.assembly_name
        yield "file_name", cfg.file_name


class UcscReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ucsc"
    source = Source.ucsc

    @dictify
    def _get_params_run(self, wildcards):
        yield from self._base_args().items()
        cfg = self.config.ucsc
        yield "db", cfg.db
        yield "file_name", cfg.file_name


class ReferenceDownloadWorkflow(BaseStep):
    name = "reference_download"
    produces = [
        DataSignature(DataType.RAW, frozenset({"ensembl", "reference", "dna"})),
        DataSignature(DataType.RAW, frozenset({"ensembl", "reference", "rna"})),
        DataSignature(DataType.RAW, frozenset({"ensembl", "reference", "protein"})),
        DataSignature(DataType.RAW, frozenset({"refseq", "reference", "dna"})),
        DataSignature(DataType.RAW, frozenset({"refseq", "reference", "rna"})),
        DataSignature(DataType.RAW, frozenset({"refseq", "reference", "protein"})),
        DataSignature(DataType.RAW, frozenset({"ucsc", "reference", "dna"})),
        DataSignature(DataType.RAW, frozenset({"ucsc", "reference", "rna"})),
        DataSignature(DataType.RAW, frozenset({"ucsc", "reference", "protein"})),
    ]

    sheet_shortcut_class = GenericSampleSheet
    config_model_class = ReferenceDownloadConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        fasta = kwargs.get("fasta", "work/reference_download/out/reference.fa")
        return {"fasta": fasta}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        source_to_class = {
            Source.ensembl: EnsemblReferenceDownloadStepPart,
            Source.refseq: RefseqReferenceDownloadStepPart,
            Source.ucsc: UcscReferenceDownloadStepPart,
        }
        self.register_sub_step_classes((source_to_class[self.config.source],))

    @listify
    def get_result_files(self):
        yield "output/reference_download/out/reference.fa"
        yield "output/reference_download/out/reference.fa.md5"

# -*- coding: utf-8 -*-
"""Implementation of the ``reference_download`` step."""

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ReferenceDownload as ReferenceDownloadConfigModel
from .model import Source

DEFAULT_CONFIG = ReferenceDownloadConfigModel.default_config_yaml_string()


class _ReferenceDownloadStepPart(BaseStepPart):
    source: Source
    actions = ("run",)

    def get_input_files(self, action):
        self._validate_action(action)
        return []

    @dictify
    def get_output_files(self, action):
        self._validate_action(action)
        output_fasta = self.config.path_output_fasta or "output/reference_download/out/reference.fa"
        work_fasta = "work/reference_download/out/reference.fa"
        yield "fasta", work_fasta
        yield "fasta_md5", work_fasta + ".md5"
        yield "output_links", [output_fasta, output_fasta + ".md5"]

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

    @dictify
    def _get_log_file(self, action):
        self._validate_action(action)
        prefix = f"work/reference_download/log/reference_download.{self.name}"
        for key, ext in (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
        ):
            yield key, prefix + ext

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield from self._base_args().items()


class EnsemblReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ensembl"
    source = Source.ensembl

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield from self._base_args().items()
        cfg = self.config.ensembl
        yield "subset", cfg.subset
        yield "file_name", cfg.file_name
        yield "datatype", cfg.datatype
        yield "chromosome", cfg.chromosome
        yield "branch", cfg.branch

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=2, runtime="2h", mem="2GB")


class RefseqReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "refseq"
    source = Source.refseq

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield from self._base_args().items()
        cfg = self.config.refseq
        yield "assembly_accession", cfg.assembly_accession
        yield "assembly_name", cfg.assembly_name
        yield "file_name", cfg.file_name

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=2, runtime="2h", mem="2GB")


class UcscReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ucsc"
    source = Source.ucsc

    @dictify
    def get_args(self, action):
        self._validate_action(action)
        yield from self._base_args().items()
        cfg = self.config.ucsc
        yield "db", cfg.db
        yield "file_name", cfg.file_name

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(threads=2, runtime="2h", mem="2GB")


class ReferenceDownloadWorkflow(BaseStep):
    name = "reference_download"
    consumes = {}
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
    def default_config_yaml(cls):
        return DEFAULT_CONFIG

    @classmethod
    def get_output_paths(cls, signature=None, **kwargs) -> dict[str, str]:
        cls.require_signature(signature)
        fasta = kwargs.get("fasta", "output/reference_download/out/reference.fa")
        return {"fasta": fasta}

    def __init__(
        self,
        workflow,
        config,
        config_lookup_paths,
        config_paths,
        workdir,
        task_name: str = "",
        **kwargs,
    ):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            task_name=task_name,
            **kwargs,
        )
        source_to_class = {
            Source.ensembl: EnsemblReferenceDownloadStepPart,
            Source.refseq: RefseqReferenceDownloadStepPart,
            Source.ucsc: UcscReferenceDownloadStepPart,
        }
        self.register_sub_step_classes((source_to_class[self.config.source],))

    def get_output_fasta_path(self) -> str:
        if self.config.path_output_fasta:
            return self.config.path_output_fasta
        return "output/reference_download/out/reference.fa"

    @listify
    def get_result_files(self):
        yield self.get_output_fasta_path()
        yield self.get_output_fasta_path() + ".md5"

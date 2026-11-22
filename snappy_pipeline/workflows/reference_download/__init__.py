# -*- coding: utf-8 -*-
"""Download of a reference sequence

The step downloads one FASTA file from Ensembl, RefSeq or UCSC.
The required ``source`` field selects the source, and the section of the same name selects the file:

- ``ensembl``: species, build, release (``version``) and ``datatype``, optionally single chromosomes (``chromosome``).
- ``refseq``: the NCBI assembly, by ``assembly_accession`` and ``assembly_name``, such as ``GCF_000001405.40`` and ``GRCh38.p14``.
  Both fields are required, so ``source: refseq`` needs the section.
- ``ucsc``: the UCSC genome database ``db``, such as ``hg38``.

``file_name`` picks another file of the release, assembly or database, and ``url`` replaces the generated URL.
``contigs`` and ``contigs_regex`` keep only some contigs of the file.
The task has one result per project, not per library.

Outputs
-------

Below ``tasks/<task>/output/``:

- ``reference_download/out/reference.fa`` with its index ``reference.fa.fai`` and sequence dictionary ``reference.dict``, each with an ``.md5`` file

Downstream tasks read the FASTA file through the key ``fasta`` of their ``reference`` field, and find the ``.fai`` and ``.dict`` next to it.
The task produces data of type ``raw`` with the tags ``reference``, its source and the ``molecule`` of the source section (``dna``, ``rna`` or ``protein``).

Examples
--------

GRCh38 from RefSeq, with the chromosomes and the mitochondrial genome (accessions ``NC_*``) only:

.. code-block:: yaml

    tasks:
      - step: reference_download
        name: genome
        config:
          source: refseq
          refseq:
            assembly_accession: GCF_000001405.40
            assembly_name: GRCh38.p14
            contigs_regex: '^NC_'
"""

from snappy_pipeline.utils import listify
from snappy_pipeline.workflows.abstract import RuleOptions, BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ReferenceDownload as ReferenceDownloadConfigModel
from .model import Source


class _ReferenceDownloadStepPart(BaseStepPart):
    actions = ("run",)
    default_resource_usage = ResourceUsage(threads=2, runtime="2h", mem="2GB")
    source: Source

    def _get_input_files_run(self, wildcards):
        return []

    def _get_output_files_run(self):
        work_fasta = "work/reference_download/out/reference.fa"
        return {"fasta": work_fasta, "fasta_md5": work_fasta + ".md5"}

    def _get_log_file_run(self):
        return {"log": f"work/reference_download/log/reference_download.{self.name}.log"}

    def _get_params_run(self, wildcards):
        # The molecule only tags the task's output; the wrapper does not read it
        return getattr(self.config, self.name).model_dump(exclude={"molecule"})


class EnsemblReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ensembl"
    wrappers = {"run": "reference_download/ensembl"}
    rule_options = {"run": RuleOptions(retries=2, cache=True)}
    source = Source.ensembl


class RefseqReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "refseq"
    wrappers = {"run": "reference_download/refseq"}
    rule_options = {"run": RuleOptions(retries=2, cache=True)}
    source = Source.refseq


class UcscReferenceDownloadStepPart(_ReferenceDownloadStepPart):
    name = "ucsc"
    wrappers = {"run": "reference_download/ucsc"}
    rule_options = {"run": RuleOptions(retries=2, cache=True)}
    source = Source.ucsc


class SidecarsStepPart(BaseStepPart):
    """The ``.fai`` and ``.dict`` next to the downloaded FASTA, which consumers declare as inputs."""

    name = "sidecars"
    actions = ("run",)
    wrappers = {"run": "reference_download/sidecars"}
    rule_options = {"run": RuleOptions(local=True)}
    default_resource_usage = ResourceUsage(threads=1, runtime="1h", mem="2GB")

    def _get_input_files_run(self, wildcards):
        return {"fasta": "work/reference_download/out/reference.fa"}

    def _get_output_files_run(self):
        prefix = "work/reference_download/out/reference"
        return {
            "fai": prefix + ".fa.fai",
            "fai_md5": prefix + ".fa.fai.md5",
            "dict": prefix + ".dict",
            "dict_md5": prefix + ".dict.md5",
        }

    def _get_log_file_run(self):
        return {"log": "work/reference_download/log/reference_download.sidecars.log"}


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

    config_model_class = ReferenceDownloadConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        molecule = getattr(config, config.source).molecule
        return (
            DataSignature(
                DataType.RAW, frozenset({str(config.source), "reference", str(molecule)})
            ),
        )

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
        self.register_sub_step_classes((source_to_class[self.config.source], SidecarsStepPart))

    @listify
    def get_result_files(self):
        yield "work/reference_download/out/reference.fa"
        yield "work/reference_download/out/reference.fa.md5"
        for name in ("reference.fa.fai", "reference.dict"):
            yield f"work/reference_download/out/{name}"
            yield f"work/reference_download/out/{name}.md5"

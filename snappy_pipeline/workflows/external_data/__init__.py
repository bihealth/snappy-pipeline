# -*- coding: utf-8 -*-
"""Existing data from outside the project

An ``external_data`` task provides existing files to other tasks.
Its ``produces`` field declares what the files are, as a data type and tags.
When the project loads, the contracts of its consumers are checked against this declaration.
The task has no rules and no ``output/`` directory.
Consumers read the files where they are.

The task has one of three modes.
With ``files``, it provides project-wide files by output key, such as a reference genome.
A key may name a list of files, such as several transcript databases of mehari.
With ``library_files``, it provides the files of each library from files tables: TSV files with one row per read group (``type: raw``: ``library``, ``readgroup``, ``left`` and ``right``) or one row per library (``library`` and the output keys).
An empty cell is no value, and relative paths are relative to the table's directory.
``mixed_se_pe: true`` accepts read groups without a ``right`` file.
``snappy scan-files`` writes the first table of each task with a ``scan`` section from the libraries of the sample sheets: the files of a library are in a folder named like its ``folderName`` in the sheet, else like the library, at any depth below one of the ``search_paths``; ``search_patterns`` map output keys to regular expressions for the path below that folder, and the ``(?P<readgroup>...)`` group of a read pattern pairs the mates.
A library without rows is an error when the project loads, and so is a row for a library that no sheet has.
With ``kits``, it provides files per enrichment kit, by kit name and output key.
A targeted library (``WES``, ``Panel_seq``) uses the kit named exactly like its ``libraryKit`` in the sample sheet.
If a library's kit is missing, its consumer fails when the project loads.
Kit names may contain letters, digits, ``_``, ``.`` and ``-``.
Every kit has the same keys.
Steps read the kit regions through ``depends_on.kits``, which needs the BED files ``targets`` and ``baits``.
If a kit has only one BED file, give it for both keys.

BED files of regions have ``type: regions``.
A single BED file has the key ``bed``, as in ``files: {bed: /refs/GRCh38/blacklist.bed}``.
Tags say what the regions are, such as ``exons``.

The output keys must be the keys that the consumer reads through its ``depends_on`` field.
Their names are checked against that field when the project loads.

Examples
--------

The reads of the libraries, a reference and a bwa index as input of ``ngs_mapping``.
``snappy scan-files`` writes ``reads.tsv`` from the ``scan`` section:

.. code-block:: yaml

    tasks:
      - step: external_data
        name: genome
        config:
          produces: {type: raw, tags: [reference, dna]}
          files: {fasta: /refs/GRCh38/GRCh38.fa}

      - step: external_data
        name: bwa_index
        config:
          produces: {type: index, tags: [bwa, dna]}
          files: {index: /refs/GRCh38/bwa/GRCh38}  # the prefix of the .amb, .ann, ... files

      - step: external_data
        name: reads
        config:
          produces: {type: raw, tags: [fastq]}
          library_files: [reads.tsv]
          scan:
            search_paths: [raw]
            search_patterns:
              - left: '(?P<readgroup>.+)\\.R1\\.fastq\\.gz'
                right: '(?P<readgroup>.+)\\.R2\\.fastq\\.gz'

      - step: ngs_mapping
        name: mapping
        config:
          depends_on:
            reads: reads
            index: bwa_index
            reference: genome
          tool: bwa
          bwa: {}

The table ``reads.tsv`` lists one read group per line (``right`` is empty for single-end data):

.. code-block:: text

    library            readgroup  left                              right
    P001-N1-DNA1-WES1  L001       raw/P001-N1-DNA1-WES1/L001.R1.fastq.gz  raw/P001-N1-DNA1-WES1/L001.R2.fastq.gz

Existing alignments as the start of variant calling, one row per library:

.. code-block:: yaml

    - step: external_data
      name: alignments
      config:
        produces: {type: alignments, tags: [dna]}
        library_files: [bams.tsv]   # the columns library, bam, bai

The regions of two exome kits, read by the coverage QC of ``ngs_data_qc``:

.. code-block:: yaml

    tasks:
      - step: external_data
        name: genome
        config:
          produces: {type: raw, tags: [reference, dna]}
          files: {fasta: /refs/GRCh38/GRCh38.fa}

      - step: external_data
        name: bwa_index
        config:
          produces: {type: index, tags: [bwa, dna]}
          files: {index: /refs/GRCh38/bwa/GRCh38}

      - step: external_data
        name: kits
        config:
          produces: {type: regions, tags: [kits]}
          kits:
            Agilent_SureSelect_V6:
              targets: /refs/GRCh38/kits/SureSelect_V6.covered.bed
              baits: /refs/GRCh38/kits/SureSelect_V6.padded.bed
            Twist_Exome_2:
              targets: /refs/GRCh38/kits/Twist_Exome_2.targets.bed
              baits: /refs/GRCh38/kits/Twist_Exome_2.probes.bed

      - name: reads
        step: external_data
        config:
          produces: {type: raw, tags: [fastq]}
          library_files: [reads.tsv]
      - step: ngs_mapping
        name: mapping
        config:
          depends_on:
            reads: reads
            index: bwa_index
            reference: genome
          tool: bwa
          bwa: {}
      - step: ngs_data_qc
        name: coverage_qc
        config:
          depends_on:
            alignments: mapping
            reference: genome
            kits: kits
          tool: alfred
          alfred: {}

Existing per-library variant calls, listed by ``snappy scan-files`` from ``/data/dragen/<library name>/``:

.. code-block:: yaml

    - step: external_data
      name: dragen_calls
      config:
        produces: {type: variants, tags: [germline, snv, indel]}
        library_files: [dragen_calls.tsv]   # the columns library, vcf, vcf_tbi
        scan:
          search_paths: [/data/dragen]
          search_patterns:
            - {vcf: '.+\\.hard-filtered\\.vcf\\.gz', vcf_tbi: '.+\\.hard-filtered\\.vcf\\.gz\\.tbi'}
"""

import os

from snappy_pipeline.base import ConfigError
from snappy_pipeline.reads import files_table_keys
from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ExternalData


def _abspath(path: str | list[str]) -> str | list[str]:
    return os.path.abspath(path) if isinstance(path, str) else [os.path.abspath(p) for p in path]


def no_rows_message(config: ExternalData, library_name: str) -> str:
    """The error for a library without rows in the files tables of ``config``, with the way to add
    them."""
    tables = ", ".join(config.library_files)
    fix = (
        "run `snappy scan-files` to add them"
        if config.scan is not None
        else "add them to the table"
    )
    return f"The files tables {tables} list no files of library {library_name!r}; {fix}"


class ExternalDataWorkflow(BaseStep):
    """Provides existing files to other tasks; has no rules of its own."""

    name = "external_data"
    produces = [DataSignature(DataType.RAW)]
    config_model_class = ExternalData

    #: The files tables are read once per run, by the project
    needs_project = True

    @classmethod
    def task_produces(cls, config, upstream):
        """The signature declared in ``produces``."""
        return (DataSignature(DataType(config.produces.type), frozenset(config.produces.tags)),)

    @classmethod
    def output_keys(cls, config, signature=None) -> set[str]:
        """The keys of ``files``, of each kit, or of the files tables (``None`` while a table
        does not exist yet)."""
        if config.files:
            return set(config.files)
        if config.kits:
            return set(next(iter(config.kits.values())))
        try:
            return files_table_keys(config.library_files, config.patterns_find_reads())
        except ValueError as e:
            raise ConfigError(str(e)) from None

    @classmethod
    def get_output_paths(
        cls,
        config,
        signature=None,
        library_name=None,
        library_kit=None,
        project=None,
        **kwargs,
    ) -> dict[str, str]:
        """Return the absolute paths of the files: project-wide, of ``library_name`` or of
        ``library_kit``."""
        if config.files:
            return {key: _abspath(path) for key, path in config.files.items()}
        if config.kits:
            if library_kit not in config.kits:
                raise ConfigError(f"No kit {library_kit!r}; known kits: {sorted(config.kits)}")
            return {key: os.path.abspath(path) for key, path in config.kits[library_kit].items()}
        files = project.files_table(config).files(library_name)
        if not files:
            raise ConfigError(no_rows_message(config, library_name))
        return files

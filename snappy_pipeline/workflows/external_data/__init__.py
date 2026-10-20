# -*- coding: utf-8 -*-
"""Existing data from outside the project

An ``external_data`` task provides existing files to other tasks.
It declares what the files are in ``produces`` (a data type and tags), so that the contracts of its consumers are checked when the project loads.
The task has no rules and no ``output/`` directory.
Consumers read the files where they are.

The task has one of three modes.
With ``files``, it provides project-wide files by output key, such as a reference genome.
A key may name a list of files, such as several transcript databases of mehari.
With ``search_paths`` and ``search_patterns``, it provides per-library files from a folder named like the library, below one of the search paths.
Each pattern maps output keys to regular expressions for the path below that folder.
Reads (``type: raw``) need the keys ``left`` (and ``right``) and a ``(?P<readgroup>...)`` group, which pairs the mates, and ``mixed_se_pe: true`` accepts read groups without a ``right`` file.
Per-library files are found like the reads of a data set, except that the folder is named like the library, not by the sample sheet.
With ``kits``, it provides files per enrichment kit, by kit name and output key.
A targeted library (``WES``, ``Panel_seq``) uses the kit named exactly like its ``libraryKit`` in the sample sheet, and a consumer fails when the project loads if a library's kit is missing.
Kit names may contain letters, digits, ``_``, ``.`` and ``-``, and every kit has the same keys.
Steps read the kit regions through ``depends_on.kits``, which needs the keys ``targets`` and ``baits`` (BED files; give one file twice if the kit has only one).

BED files of regions have ``type: regions``.
A single BED file has the key ``bed``, as in ``files: {bed: /refs/GRCh38/blacklist.bed}``, and tags say what the regions are, such as ``exons``.

The output keys must be the keys that the consumer reads through its ``depends_on`` field.
Their names are checked against that field when the project loads.

Examples
--------

Existing reads, a reference and a bwa index as input of ``ngs_mapping``:

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
        name: trimmed_reads
        config:
          produces: {type: raw, tags: [trimmed]}
          search_paths: [/data/earlier_run/trimmed]
          search_patterns:
            - left: '(?P<readgroup>.+)_R1\\.fastq\\.gz'
              right: '(?P<readgroup>.+)_R2\\.fastq\\.gz'

      - step: ngs_mapping
        name: mapping
        config:
          depends_on:
            reads: trimmed_reads
            index: bwa_index
            reference: genome
          tool: bwa
          bwa: {}

The regions of two exome kits, read by the coverage report of ``ngs_mapping``:

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

      - step: ngs_mapping
        name: mapping
        config:
          depends_on:
            reads: data_sets
            index: bwa_index
            reference: genome
            kits: kits
          tool: bwa
          bwa: {}
          target_coverage_report:
            enabled: true

Existing per-library variant calls, found as ``<search path>/<library name>/...``:

.. code-block:: yaml

    - step: external_data
      name: dragen_calls
      config:
        produces: {type: variants, tags: [germline, snv, indel]}
        search_paths: [/data/dragen]
        search_patterns:
          - {vcf: '.+\\.hard-filtered\\.vcf\\.gz', vcf_tbi: '.+\\.hard-filtered\\.vcf\\.gz\\.tbi'}
"""

import os

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ExternalData


def _abspath(path: str | list[str]) -> str | list[str]:
    return os.path.abspath(path) if isinstance(path, str) else [os.path.abspath(p) for p in path]


class ExternalDataWorkflow(BaseStep):
    """Provides existing files to other tasks; has no rules of its own."""

    name = "external_data"
    produces = [DataSignature(DataType.RAW)]
    sheet_shortcut_class = GenericSampleSheet
    config_model_class = ExternalData

    #: Per-library files are found with the project's read discovery
    needs_discovery = True

    @classmethod
    def task_produces(cls, config, upstream):
        """The signature declared in ``produces``."""
        return (DataSignature(DataType(config.produces.type), frozenset(config.produces.tags)),)

    @classmethod
    def output_keys(cls, config, signature=None) -> set[str]:
        """The keys of ``files``, of the search patterns or of each kit."""
        if config.files:
            return set(config.files)
        if config.kits:
            return set(next(iter(config.kits.values())))
        return {key for pattern in config.search_patterns for key in pattern}

    @classmethod
    def get_output_paths(
        cls, config, signature=None, library_name=None, library_kit=None, discovery=None, **kwargs
    ) -> dict[str, str]:
        """Return the absolute paths of the files: project-wide, of ``library_name`` or of
        ``library_kit``."""
        if config.files:
            return {key: _abspath(path) for key, path in config.files.items()}
        if config.kits:
            if library_kit not in config.kits:
                raise ValueError(f"No kit {library_kit!r}; known kits: {sorted(config.kits)}")
            return {key: os.path.abspath(path) for key, path in config.kits[library_kit].items()}
        roots = [os.path.abspath(path) for path in config.search_paths]
        return discovery.find_files(roots, library_name, config.search_patterns)

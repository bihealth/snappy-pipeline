# -*- coding: utf-8 -*-
"""Implementation of the ``external_data`` step

An ``external_data`` task provides existing data from outside the project to other tasks. It
declares what the data is (``produces``), so contracts check its consumers when the project
loads, and it has no rules of its own: consumers read the files where they are.

Two modes:

- ``files``: project-wide files by output key, e.g. a reference with its sidecars.
- ``search_paths`` and ``search_patterns``: per-library data, in a folder named like the library
  below one of the search paths. Each pattern maps output keys to regular expressions for the
  path below that folder. Reads (``type: raw``) need ``left`` (and ``right``) keys and a
  ``(?P<readgroup>...)`` group, which pairs the mates.

Example task config
-------------------

.. code-block:: yaml

    tasks:
      - step: external_data
        name: dragen_calls
        config:
          produces: {type: variants, tags: [germline, snv, indel]}
          search_paths: [/data/dragen]
          search_patterns:
            - {vcf: '.+\\.hard-filtered\\.vcf\\.gz', vcf_tbi: '.+\\.hard-filtered\\.vcf\\.gz\\.tbi'}

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
          tool: bwa
          bwa:
            path_index: /path/to/bwa/index.fa
"""

import os

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.workflows.abstract import BaseStep
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import ExternalData


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
    def output_keys(cls, config) -> set[str]:
        """The keys of ``files``, or the keys of the search patterns."""
        if config.files:
            return set(config.files)
        return {key for pattern in config.search_patterns for key in pattern}

    @classmethod
    def get_output_paths(
        cls, config, signature=None, library_name=None, discovery=None, **kwargs
    ) -> dict[str, str]:
        """Return the absolute paths of the files, project-wide or of ``library_name``."""
        if config.files:
            return {key: os.path.abspath(path) for key, path in config.files.items()}
        roots = [os.path.abspath(path) for path in config.search_paths]
        return discovery.find_files(roots, library_name, config.search_patterns)

    def get_result_files(self):
        """No output files: the data already exists."""
        return []

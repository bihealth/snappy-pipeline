# -*- coding: utf-8 -*-
"""Implementation of the ``link_in`` step

The ``link_in`` step is a lightweight configuration-carrier task that declares an external
directory containing pre-processed FASTQ files to be linked into a downstream workflow.

It produces no output files of its own. Its sole purpose is to publish a ``path`` that
other workflows (e.g. ``ngs_mapping``, ``adapter_trimming``, ``hla_typing``) can resolve
via ``depends_on: {link_in: <task_name>}`` to override the default ``data_sets``-based
FASTQ search path.

Example task config
-------------------

.. code-block:: yaml

    tasks:
      - step: link_in
        name: trimmed_reads
        config:
          path: /path/to/adapter_trimming/output/bbduk

      - step: ngs_mapping
        name: ngs_mapping
        depends_on:
          link_in: trimmed_reads
        config:
          tool: bwa
          bwa:
            path_index: /path/to/bwa/index.fa

"""

from biomedsheets.shortcuts import GenericSampleSheet

from snappy_pipeline.workflows.abstract import BaseStep

from .model import LinkIn


class LinkInWorkflow(BaseStep):
    """Config-carrier step declaring an external preprocessed FASTQ source directory.

    Use this as an upstream dependency for any workflow step that needs to link in
    pre-processed FASTQs from outside the current pipeline run, instead of crawling
    the ``data_sets`` search paths.
    """

    name = "link_in"
    sheet_shortcut_class = GenericSampleSheet
    config_model_class = LinkIn

    def get_result_files(self):
        """No output files — this step only carries configuration."""
        return []

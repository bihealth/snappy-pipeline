.. _installation:

============
Installation
============

.. note::

    If you are on the BIH cluster, first read :ref:`quickstart` as this also explains the temporary directory.

-------------
Prerequisites
-------------

Install `pixi <https://pixi.sh>`_ (see https://pixi.sh/latest/#installation).
The CUBI pipeline uses pixi to manage all dependencies -- both conda packages
(system tools like BWA, STAR, samtools) and PyPI packages.

For cluster execution, you need a Snakemake profile available (use ``--slurm``
with the ``snappy run`` command if your cluster uses SLURM).

-------------------------
User Installation
-------------------------

If you just want to *run* a pipeline (not develop it), clone the repository and
let pixi create the environment:

.. code-block:: shell

    $ git clone git@github.com:bihealth/snappy-pipeline.git
    $ cd snappy-pipeline
    $ pixi install

After installation the ``snappy`` command is available via
``pixi run snappy <subcommand> ...``, or by activating the environment with
``eval "$(pixi shell-hook)"``.

For a reproducible install pinned to the exact dependency versions in the
lock file, use ``pixi install --frozen`` instead.

-------------------------
Installing as a Developer
-------------------------

Same clone + install steps as above, then use the ``dev`` environment which
includes test, lint, and documentation tools:

.. code-block:: shell

    $ pixi install
    $ pixi run -e dev --  # one-off commands
    $ pixi shell -e dev   # activate dev environment

Running the Tests
=================

.. code-block:: shell

    $ pixi run -e dev test

Running the Style Checks
=========================

.. code-block:: shell

    $ pixi run -e dev lint             # ruff check + ruff format --check + snakefmt --check
    $ pixi run -e dev fmt              # auto-format with ruff
    $ pixi run -e dev snakefmt         # auto-format Snakemake files
    $ pixi run -e dev srcfmt           # run all formatters

Building the Documentation
==========================

.. code-block:: shell

    $ pixi run -e docs docs

The HTML output is written to ``docs/_build/html/``.

Developer Documentation
=======================

Make sure to also read the "Pipeline Developer Docs" section, starting with :ref:`dev_intro`.

Configuring GATK3
==================

Some wrappers rely on GATK 3.
GATK v3 is not free software and cannot be redistributed.
If you are a member of CUBI, you can use the central GATK download.
Alternatively, download the tarball `from the Broad archive <https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2>`_.

To register GATKv3 with the pipeline, create the conda environments first, then
register the tarball into each environment that requires it:

.. code-block:: shell

    $ cd /path/to/project
    $ snappy run -- --conda-create-envs-only
    $ grep 'gatk.*3' .snakemake/conda/*.yaml
    .snakemake/conda/d76b719b718c942f8e49e55059e956a6.yaml:  - gatk =3
    $ for yaml in $(grep -l 'gatk.*3' .snakemake/conda/*.yaml); do
          environ=${yaml%.yaml}
          conda activate $environ
          gatk3-register /path/to/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
          conda deactivate
      done

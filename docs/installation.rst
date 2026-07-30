.. _installation:

============
Installation
============

.. note::

    If you are on the BIH cluster, first read :ref:`quickstart` as this also explains the temporary directory.

-------------
Prerequisites
-------------

The CUBI pipeline requires Python >=3.12 (e.g., from a pixi or Miniconda3 installation).

For cluster execution, you need a Snakemake profile available (use ``--slurm`` with the ``snappy run`` command if your cluster uses SLURM).

-------------------------
Installing as a Developer
-------------------------

We use `pixi <https://pixi.sh>`_ as the project manager (install pixi: see https://pixi.sh/latest/#installation).
Pixi handles both conda dependencies (system tools like BWA, STAR, samtools) and PyPI dependencies.

.. code-block:: shell

    $ git clone git@github.com:bihealth/snappy-pipeline.git
    $ cd snappy-pipeline
    $ pixi install

This sets up all environments, including the ``dev`` environment with test and linting tools.

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

Developer Documentation
=======================

Make sure to also read the "Pipeline Developer Docs" section, starting with :ref:`dev_intro`.

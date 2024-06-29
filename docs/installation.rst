.. _installation:

============
Installation
============

.. note::

    If you are on the BIH cluster, first read :ref:`quickstart` as this also explains the temporary directory.

-------------
Prerequisites
-------------

The CUBI pipeline requires Python >=3.12 (e.g., from a Miniconda3 installation).

More recent versions also work but other requirements as Snakemake might make it depend on a more recent Python version.

For cluster execution, you need a Snakemake profile available.

--------------------
Installing a Release
--------------------

This is the recommended way if you just want to use the pipeline, simply read :ref:`quickstart`.

-------------------------
Installing as a Developer
-------------------------

It is highly recommended to have a Miniconda installation for the development as this allows for easily resetting everything.
You can of course clone the code anywhere you like.

.. code-block:: shell

    $ mkdir -p ~/Development/pipeline_dev
    $ cd ~/Development/pipeline_dev
    $ git clone git@github.com:bihealth/snappy-pipeline.git
    $ cd snappy_pipeline
    $ conda env create -n snappy_dev --file environment.yaml
    $ conda activate snappy_dev
    $ pip install -e ".[all]"


Installing pre-commit-hooks
===========================
To make it easier to follow the coding style, we use `pre-commit <https://pre-commit.com>`_ hooks.
These hooks will run the style checks before you commit your changes and will automatically fix some issues.

First, install the pre-commit package (if not already installed):

.. code-block:: shell

    $ conda install pre-commit  # or pip install pre-commit

Then, install the pre-commit hooks:

.. code-block:: shell

    $ pre-commit install

The next time you commit changes, the pre-commit hooks will run automatically.

Running the Tests
=================

To run the tests, simply invoke ``pytest``:

.. code-block:: shell

    $ cd ~/Development/pipeline_dev
    $ pytest

Running the Style Checks
========================

.. code-block:: shell

    $ cd ~/Development/pipeline_dev
    $ make lint

Developer Documentation
=======================

Make sure to also read the "Pipeline Developer Docs" section, starting with :ref:`dev_intro`.

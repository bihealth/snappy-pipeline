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
    $ pip install -e .
    $ pip install -r requirements/dev.txt

It's also a good idea to install some packages required for testing through conda:

.. code-block:: shell

    $ conda env update --name root --file environment.yaml

(If you do not do this, please make sure that you have git-lfs in your PATH through other means)

Running the Tests
=================

To run the tests, you need to add the packages in ``requirements/test.txt``.

.. code-block:: shell

    $ cd ~/Development/pipeline_dev
    $ py.test

Running the Style Checks
========================

.. code-block:: shell

    $ cd ~/Development/pipeline_dev
    $ flake8

Developer Documentation
=======================

Make sure to also read the "Pipeline Developer Docs" section, starting with :ref:`dev_intro`.

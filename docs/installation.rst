.. _installation:

============
Installation
============

.. note::

    If you are on the BIH cluster, first read :ref:`quickstart` as this also explains how to setup DRMAA and the temporary directory.

-------------
Prerequisites
-------------

The CUBI pipeline requires Python >=3.4 (e.g., from a Miniconda3 installation).

More recent versions also work but other requirements as Snakemake might make it depend on a more recent Python version.

For cluster execution, you need:

- a working DRMAA installation and proper configuration of ``DRMAA_LIBRARY_PATH`` environment variable
- the Python package ``drmaa``

--------------------
Installing a Release
--------------------

This is the recommended way if you just want to use the pipeline, simply read :ref:`quickstart`.

-------------------------
Installing as a Developer
-------------------------

It's quite straightforward, just make sure that you have properly setup the ``DRMAA_LIBRARY_PATH`` environment variable when you're on the cluster.
Also, it is highly recommended to have a Miniconda installation for the development as this allows for easily resetting everything.
You can of course clone the code anywhere you like.

.. code-block:: shell

    $ mkdir -p ~/Development/pipeline_dev
    $ cd ~/Development/pipeline_dev
    $ git clone git@gitlab.bihealth.org:cubi/snappy_pipeline.git
    $ cd snappy_pipeline
    $ pip install -e .
    $ pip install -r requirements/dev.txt

It's also a good idea to install some packages required for testing through conda:

.. code-block:: shell

    $ conda env update --name root --file environment.yaml

(If you do not do this, please make sure that you have git-lfs in your PATH through other means.)

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

.. _quickstart:

==========
Quickstart
==========

This chapter gives the minimal number of commands required for setting up the pipeline on the BIH cluster.

.. note::

    This describes the setup as a pipeline user.
    If you want to know about the setup as a pipeline developer, see :ref:`installation`.

-----------------------
Install pixi
-----------------------

First, install `pixi <https://pixi.sh>`_ (see https://pixi.sh/latest/#installation for instructions).

-----------------------
Install Snappy Pipeline
-----------------------

The recommended way of installing the CUBI pipeline is via pip inside a pixi environment.

.. code-block:: shell

    $ VERSION=vX.Y.Z
    $ pip install git+ssh://git@github.com:bihealth/snappy-pipeline.git@v${VERSION}#egg=snappy_pipeline

Or see ``README.rst`` for a more detailed installation guide.

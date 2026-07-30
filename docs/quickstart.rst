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

Clone the repository and install with pixi:

.. code-block:: shell

    $ git clone git@github.com:bihealth/snappy-pipeline.git
    $ cd snappy-pipeline
    $ pixi install

To pin to the exact versions in the lock file (``pixi.lock``), use
``pixi install --frozen``. The lock file is the single source of truth for all
dependency versions and is kept in sync with ``pyproject.toml`` via
``pixi update``.

After installation the ``snappy`` command is available via
``pixi run snappy <subcommand> ...``.

Or see :ref:`installation` for a more detailed guide, including the developer
setup with test, lint, and documentation tooling.

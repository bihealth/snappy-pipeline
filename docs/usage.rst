.. _usage:

=====
Usage
=====

As a user, you will mostly interface with the CUBI pipeline system using the ``snappy-snake`` program.

This program is a wrapper around `Snakemake <https://snakemake.bitbucket.org>`_ and provides the following features:

- a number of pre-packaged, well-tested workflows (pipeline steps) that are
- driven by configuration and sample sheet files and
- can be shared over multiple projects; it has
- good integration into DRMAA-based scheduling systems and a
- easy-to-use command line interface.

Here is how to get command line help:

.. code-block:: shell

    $ snappy-snake --help

Here is how to print the version and enabled features:

.. code-block:: shell

    $ snappy-snake --snappy-pipeline-self-test

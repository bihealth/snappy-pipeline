.. _usage:

=====
Usage
=====

As a user, you will mostly interface with the CUBI pipeline system using the ``snappy`` command.

All subcommands hang off the ``snappy`` root group:

.. code-block:: shell

    $ snappy --help

Subcommands
===========

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Subcommand
     - Purpose
   * - ``snappy init``
     - Scaffold a new project -- creates ``config.yaml``, ``pipeline_job.sh``, and a sample sheet template.
   * - ``snappy task add``
     - Add a new task to an existing project's ``config.yaml``.
   * - ``snappy task list``
     - List all tasks defined in the project.
   * - ``snappy run``
     - Run the pipeline via Snakemake.  This is the main entry point for executing workflows.
   * - ``snappy watch``
     - Launch the ``snkmt`` TUI to monitor a running workflow via its SQLite database.
   * - ``snappy refresh``
     - Recreate ``pipeline_job.sh`` and ensure ``slurm_log`` exists.
   * - ``snappy status``
     - Check the status of a SLURM job via ``sacct``.
   * - ``snappy pull-sheet``
     - Pull a sample sheet from the SODAR API.

snappy run options
==================

``--task <name>``
    Build only the named task (overrides the default leaf-task selection).

``--all-tasks``
    Build every task in the config.

``--slurm``
    Enable SLURM cluster execution (layers a SLURM profile on top of the
    default conda profile).

Everything after ``--`` is passed verbatim to Snakemake.

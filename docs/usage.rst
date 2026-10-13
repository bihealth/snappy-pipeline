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
   * - ``snappy report``
     - Build the Snakemake report of the QC outputs (``report.zip`` in the project directory).
   * - ``snappy logs``
     - Archive the logs of a run whose Snakemake process was killed (``snappy run`` does this
       itself when it ends).
   * - ``snappy refresh``
     - Recreate ``pipeline_job.sh`` and ensure ``slurm_log`` exists.
   * - ``snappy status``
     - Check the status of a SLURM job via ``sacct``.
   * - ``snappy pull-sheet``
     - Pull a sample sheet from the SODAR API.

snappy run options
==================

``--task <name>``
    Build only the named task (overrides the default leaf-task selection). Upstream tasks
    whose outputs are missing or outdated are built as well.

``--frozen``
    With ``--task``: load only that task's rules. Outputs of upstream tasks are used as they
    are and never rebuilt; if one is missing, Snakemake stops with a missing-input error that
    names the file.

``--all-tasks``
    Build every task in the config.

``--slurm``
    Enable SLURM cluster execution (layers a SLURM profile on top of the
    default conda profile).

``-v``, ``--verbose``
    Increase verbosity; prints the resolved task configuration and the
    Snakemake command line.

Everything after ``--`` is passed verbatim to Snakemake.

Report
======

``snappy report`` builds Snakemake's report of the project: the QC outputs of the tasks, one
category per task (FastQC, Picard metrics, ``bcftools stats``, cnvkit plots). It uses the same
Snakefile, directory and profile as ``snappy run``; arguments after ``--`` go to Snakemake. The
default ``report.zip`` can hold HTML directories such as the FastQC reports.

Log archive
===========

When a run ends, successfully or not, ``snappy run`` writes ``logs/<run start>.tar.gz``. It holds
the run's Snakemake log, the logs of the jobs that ran (with their SLURM logs), and a copy of the
``snkmt`` database. The run's Snakemake log names these files, so the task directories are not
searched. Dry-runs write no archive. If the Snakemake process was killed (walltime,
``scancel``), its end hooks did not run; ``snappy logs`` builds the archive from the newest
Snakemake log, or from the one given as argument.

Temporary files
===============

Intermediate files that only later rules of the same task read, such as the per-region outputs
of a scattered caller, are ``temp()`` files: Snakemake deletes them once every job that reads
them is done. Files in ``output/``, the files a task passes to other tasks, and project-wide
resources (installed software, prepared references) are never deleted, so ``--task`` and
``--frozen`` runs always find what they need. ``gatk4_hc_gvcf`` keeps its gVCFs for joint
re-genotyping unless ``keep_gvcfs: false``.

Once a temporary file is gone, a change to a later rule (a parameter, its code) also reruns the
rules that produced the file, for example all scatter jobs of a caller.

Snakemake options for temporary files (after ``--``):

``--notemp``
    Keep all temporary files of this run, for debugging.

``--all-temp``
    Treat every output as temporary, to save space in test runs.

``--delete-temp-output``
    Delete temporary files left over by earlier runs; with ``--dry-run``, only list them.

``--keep-incomplete``
    Keep the outputs of failed jobs, for debugging.

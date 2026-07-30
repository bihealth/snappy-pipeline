.. _overview:

========
Overview
========

This chapter gives you the big picture of the CUBI pipeline system.
The audience is people who already have experience with Bioinformatics pipeline/workflow systems and see the benefit of such systems (e.g., GNU Make, Snakemake, bpipe, etc.) over shell files over interactive bash commands.
You are part of the audience if you agree that automation is key for effective, efficient, and reproducible Bioinformatics analysis as this is a requirement for important key requirements such as provenance tracking.

Up to a certain point, automation in Bioinformatics workflows is a no-brainer as the same steps always repeat themselves.
After this point, the tasks might become very project specific and not benefit from generic, shared automation much.
One example is report generation where most of the code cannot be re-used in different projects.
Here, different means should be used (e.g., using Rmarkdown documents).

The CUBI pipeline system is aimed at the steps upstream of this "certain point".


.. _motivation:

----------
Motivation
----------

Generally, the aim was to achieve the following properties in a pipeline system:

Re-use.
    Ability to re-use common Bioinformatics analysis steps.
    Mostly, these are shell snippets with calls to standard Bioinformatics tools with some glue and conversion code thrown in.

Configurability.
    Allow for good configuration by configuration files.
    No paths should be hard-coded in the system but instead come from a configuration file.
    Further, the important parameters that might need tweaking should be exposed through the configuration.

Sensible Default Parameters.
    Provide sensible defaults for configuration.
    Ideally, use auto-tuning of parameters (e.g., call BWA-ALN for short and single reads, BWA-MEM for long, paired reads).

Good Documentation.
    Provide good documentation of the pipeline system.
    Widespread re-use improves the pay-off of good documentation.

Logging software versions.
    Log the version of the pipeline and tools to allow analyses to be repeated with the same program versions in the future.
    At the very least, knowing the versions used can help explain (slight) differences in results.

Versioning of pipeline code.
    Use semantic versioning for result files.
    Output paths should not change or disappear between minor versions.

Robustness.
    Pipeline execution failure should be prevented (e.g., all required parameters to called tools should be present) and technical weaknesses should be worked around (e.g., by allowing restarting of jobs).

Restartability.
    If the pipeline is stopped or when new input data sets are added, do not repeat unnecessary work.
    Further, if an intermediate file changes, the dependent files should be updated.
    (This is similar to what GNU Make does.)

Ease of use.
    Help the users not shoot themselves in the foot too badly (e.g., prevent accidentally overwriting already existing files).
    Easy local and cluster execution.
    At least provide sensible defaults for resource requirements, ideally auto-configured from input data.


.. _definitions:

-----------
Definitions
-----------

For clarity, this documentation uses the following definitions for separating the code for pipeline steps and the actual execution of code.

pipeline
    Code for performing a set of Bioinformatics tasks in an automated fashion.

project
    A project corresponds to a directory in the file system.
    A project is an **instance** of a pipeline, in that the different available pipeline parts are plugged together by configuration and the executed.

(pipeline) step
    Program code (Snakefiles, scripts etc.) for performing a certain "encapsulated" set of tasks.
    Examples are read mapping, variant calling, and variant annotation.

(pipeline) step instance
    A project's folder on the file system, with *configuration*, where a pipeline step is executed.
    The instance shares the pipeline step code with all other intances of the same type.

working directory
    A directory on the disk for a step instance.


.. _pipeline_projects:

------------------------
A Simple Example Project
------------------------

The above part of this chapter is quite abstract.
Let us draw some pictures and go from the abstract description to a concrete example.
We will use a simple NGS somatic variant calling pipeline for matched tumor/normal pairs, setup for WES or WGS processing.

Components of a CUBI Pipeline Project
=====================================

The following figure shows the different components that are involved for running the CUBI pipeline.

.. figure:: figures/overview_locations.*

    Overview of the different components for running the CUBI pipeline.
    Boxes of the same color indicate that the represented entities belong together.

The different parts are as follows

- The blue-colored boxes represent the ``snappy_pipeline`` Python package that contains the ``snappy`` command and the code for the different pipeline steps.

- The yellow-colored boxes represent the project directory with the ``tasks/`` subdirectory.
  Each task is namespaced below ``tasks/<task_name>/``.

- The orange-colored boxes represent the configuration.
  There is a project-wide ``config.yaml`` file that defines all tasks under a ``tasks`` key.
  Each task carries its own self-contained configuration, including the tool, library selection, and dependencies.

- The purple-colored box represents static data such as the reference sequence, annotations, databases such as dbSNP or dbNSFP.
  These static data files are created and maintained independently of the individual projects.

- The green box represents the raw input data, e.g., a directory containing the FASTQ reads for each sample.
  While, of course, raw data can be shared over projects, the data directories are usually under control of the project manager while the static data is under control of the maintainer of the static data project of **Cubit**.

- The brown box represents the bio-medical sample sheets with metadata that describe the data sets of the experiment and also (at least) parts of the experimental setup.

The number of steps might seem intimidating at first, but you will quickly get used to this arrangement.
After all, the configuration is closely related to the directories.
Further, static data and raw data paths are just put into the configuration once and otherwise you do not have to deal with it.
Also, there is UI support for generating and updating the bio-medical sample sheet files.

Components of a Pipeline Step Instance Excecution
=================================================

The following figure shows the components involved when executing a pipeline step (in this case, the NGS read mapping step).

.. figure:: figures/components_step_instance.*

    Overview of the components involved when executing a pipeline step in a working directory.

The different parts are as follows:

- The working directory ``project/ngs_mapping``.
- The step-level configuration in ``project/ngs_mapping/config.yaml``.
- The project-level configurations in ``project/.snappy_pipeline/config.yaml`` (by convention).
- The ``snappy_pipeline`` Python package installed centrally.
- The bio-medical sample sheets with the data sets to use.
  (The project-wide configuration files point at these files.)
- The static data files setup by the Cubit administrator (here, it would be the reference FASTA path and the read mapper index location).
- The raw data files to be processed by the pipeline step (here, it would be the sample FASTQ files).

How FASTQ files are found
=========================

In its ``data_sets`` section, the project-level configuration file provides search paths and search patterns to find the input FASTQ files. ``snappy`` internally combines these paths & search patterns with the sample-specific path information provided in the sample sheet. In the end, FASTQ files retained for processing are files which paths match:

::

    <configuration search path>/<sample-specific folder>/../<search pattern>

The search will loop over provided search paths & search patterns. Paired reads files are coupled by similarity of their path. Note that when the ``Folder`` entry is absent from the sample sheet, the library name is used instead.


Overview of the Somatic Variant Pipeline
========================================

The following figure shows an overview the simple somatic variant calling pipeline used in the example.

.. figure:: figures/overview_somatic_varcall.*

    Overview of the steps in somatic variant calling pipeline.

The configuration, static data files, and bio-medical sample sheets are used for the input of all pipeline steps.
The raw data files are used for the input of the NGS mapping.
The resulting read alignments are used as the input for the somatic variant calling.
The resulting somatic variant files are then used as the input for the somatic variant annotation.

Within each step the following actions are performed:

1. The reads are first mapped to a reference genome, yielding BAM files contaning the read alignments. (Additional text files with the alignment reports are also generated at this step, but this pipeline does not use these files in the downstream steps.)
2. Then, the pairs of BAM alignments for the matched tumor/normal samples for each individual are given to a somatic variant caller that produces a VCF file with the list of somatic variants for each patient.
3. Finally, variant annotations are added to indicate whether each event is present in the snp databases specified in the configuration (e.g., dbSNP or COSMIC) and functional mutation impact predictions are also added using the tool specified in the configuration (e.g., using MutationTaster).

The Matched Cancer Data Schema
==============================

For the somatic variant calling, the matched cancer study bio-medical data sheet schema is used.
It is described in full in the BioMed Sheets project.
Here, we give a summary so this document is self-contained.

- The study contains a number of patients/donors, and each individual is associated with a normal and a tumor sample.
- From each sample, an WES library is generated and sequenced; for each library, there is a directory with the library name, storing the FASTQ files from sequencing.

Project Directory Setup
=======================

A project directory is set up with the ``snappy init`` command:

.. code-block:: shell

    $ snappy init --directory somatic_project
    $ tree -a somatic_project
    somatic_project/
    +-- config.yaml
    +-- pipeline_job.sh
    +-- samplesheet.tsv
    +-- raw/
    +-- resources/

The ``config.yaml`` file contains a ``tasks`` list with the step configurations, ``static_data_config`` for paths to reference data, and a ``data_sets`` section describing the input data and sample sheet.

Tasks reference step types (``ngs_mapping``, ``variant_calling``, etc.) and each task can have its own tool, library selection, and dependency configuration.
For example, after editing ``config.yaml`` to add a mapping and variant calling task, the project might look like:

.. code-block:: yaml

    static_data_config:
      reference:
        path: ../../resources/refs/GRCh38.fa

    tasks:
      - name: bwa_mapping
        step: ngs_mapping
        config:
          tool: bwa
          bwa:
            path_index: ../../resources/refs/bwa_index

      - name: strelka_calling
        step: variant_calling
        config:
          tool: strelka
          depends_on:
            ngs_mapping: bwa_mapping

    data_sets:
      batch1:
        file: samplesheet.tsv
        search_patterns:
          - { left: '*.R1.fastq.gz', right: '*.R2.fastq.gz' }
        search_paths:
          - raw
        type: matched_cancer

Adding tasks can also be done incrementally with ``snappy task add``:

.. code-block:: shell

    $ snappy task add --directory somatic_project ngs_mapping=bwa_mapping
    $ snappy task add --directory somatic_project variant_calling=strelka_calling

Path Resolution for Static Data and Configuration Files
========================================================

Relative paths in ``static_data_config`` and in task-level configuration
are resolved relative to the **config file's directory** at validation time.
Absolute paths are supported and will not be modified.

Working Directory Layout
========================

After running a task, the output and working files are namespaced below ``tasks/<task_name>/``:

::

    somatic_project/
    +-- config.yaml
    +-- tasks/
    |   +-- bwa_mapping/
    |   |   +-- output/
    |   |   +-- work/
    |   +-- strelka_calling/
    |       +-- output/
    |       +-- work/
    +-- raw/
    +-- resources/

Adding Sample Sheets
====================

Sample sheets are TSV files describing the study design.
For matched cancer studies, the format is:

.. code-block:: tsv

    [Metadata]
    schema          cancer_matched
    schema_version  v1

    [Data]
    patientName sampleName  isTumor    libraryType folderName
    P001    N1  N   WES P001-N1-DNA1-WES1
    P001    T1  Y   WES P001-T1-DNA1-WES1
    P002    N1  N   WES P002-N1-DNA1-WES1
    P002    T1  Y   WES P002-T1-DNA1-WES1

Executing the Pipeline
======================

To run the pipeline, use ``snappy run``:

.. code-block:: shell

    $ cd somatic_project
    $ snappy run                    # run all leaf tasks
    $ snappy run --task bwa_mapping  # run only bwa_mapping

For cluster execution with SLURM:

.. code-block:: shell

    $ snappy run --slurm

Extra Snakemake arguments can be passed after ``--``:

.. code-block:: shell

    $ snappy run --task bwa_mapping -n -- --quiet

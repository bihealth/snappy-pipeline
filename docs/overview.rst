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

------------------
An Example Project
------------------

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

- The blue-colored boxes represent the ``snappy_pipeline`` Python package that contains the ``cubi-snake`` executable and the code for the different pipeline steps.

- The yellow-colored boxes represent the project directory with the different sub directories for the step instances.
  For each step that is to be executed (with a given configuration set), a directory is created.
  In the given example, there is only one directory (and thus instance) for each step.

- The orange-colored boxes represent the configuration.
  There is a project-wide ``config.yaml`` file that defines project-wide defaults.
  Each step instance can then override certain settings, similar to how sub-classing in OOP works.
  One read mapping step instance may use GRCh37 for the reference and another instance might use GRCh38 (not shown in this example).

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

The project directory is setup with the following helper tool:

.. VS Code highlighting broken with backtick in code-block, thus double-colon

::

    $ snappy-start-project --directory somatic_project
    [...]
    Do not forget to fill out your README.md file!

    SUCCESS: all done, have a nice day!

    $ tree -a somatic_project
    somatic_project/
    +-- .snappy_pipeline/
    |   `-- config.yaml
    `-- README.md

The ``config.yaml`` file is setup with common configuration for the pipeline steps.
The template used uses the paths specific to the Cubit installation on the BIH cluster.
In the far future, custom templates will be used for this and the generic files will contain "TODO" entries for changes.

Further, a project-wide ``README.md`` file is setup in which you can place documentation on the project.

.. code-block:: shell

    $ cd somatic_project
    $ head .snappy_pipeline/config.yaml
    # CUBI Pipeline Project "somatic_project" Configuration
    #
    # created: 2017-02-03T12:57:17.302044

    # Step Configuration ==============================================================================
    #
    # Configuration for paths with static data.  This has been preconfigured for the paths on the BIH
    # cluster.
    #
    static_data_config:

Working Directories for Step Instances
======================================

Next, we create the different step instances that we want to use using ``snappy-start-step``.
Note that this will extend the ``.snappy_pipeline/config.yaml`` file if there is no configuration entry for the given step.
A different name for the instance can be given using the ``--step`` parameter.

Adding the ``ngs_mapping`` step creates the required directory and configuration files pointing to the global configuration for extension.
Note how the difference in the project-wide configuration (and all other files created or modified) is displayed in the script's output.

See :ref:`step_ngs_mapping` for the default configuration of the ``ngs_mapping`` step.
For all configuration settings that have no default and are marked with a ``# required`` comment (case insensitive), these markers are copied to the project configuration so you know which settings to adjust.

::

    $ cd somatic_project
    $ snappy-start-step --step ngs_mapping
    [...]
    INFO: applying the following change:

    --- a/.snappy_pipeline/config.yaml	2017-02-03T12:47:32.246833
    +++ b/.snappy_pipeline/config.yaml	2017-02-03T12:49:29.811706
    @@ -22,7 +22,12 @@
     # Configuration for the individual steps.  These can be filled by the snappy-start-step command
     # or initialized already with snappy-start-project.
     #
    -step_config: {}
    +step_config:
    +  ngs_mapping:
    +    bwa:
    +      path_index:  # REQUIRED
    +    star:
    +      path_index:  # REQUIRED

     # Data Sets =======================================================================================
     #
    [...]

    $ tree ngs_mapping
    ngs_mapping/
    |-- config.yaml
    |-- pipeline_job.sh
    `-- sge_log

    $ cat ngs_mapping/config.yaml
    pipeline_step:
    name: ngs_mapping
    version: 1

    $ref: 'file://../.snappy/config.yaml'

Similarly, adding ``somatic_variant_calling`` adds configuration for somatic variant calling.

::

    $ snappy-start-step --step somatic_variant_calling
    [...]
    INFO: applying the following change:

    --- a/.snappy/config.yaml	2017-02-03T13:11:10.023648
    +++ b/.snappy/config.yaml	2017-02-03T13:11:20.806588
    @@ -29,6 +29,10 @@
         star:
         path_index: REQUIRED  # REQUIRED

    +  somatic_variant_calling:
    +    path_ngs_mapping: ../ngs_mapping  # REQUIRED
    +    scalpel:
    +      path_target_regions:  # REQUIRED
     # Data Sets =======================================================================================
     #
     # Define data sets.  The search paths and patterns are given per data set.
    [...]

    $ tree somatic_variant_calling
    somatic_variant_calling
    +-- sge_log/
    `-- config.yaml

The same is true for adding ``somatic_variant_annotation``.

::

    $ snappy-start-step --step somatic_variant_annotation
    [...]
    INFO: applying the following change:

    --- a/.snappy_pipeline/config.yaml	2017-02-03T13:11:20.807090
    +++ b/.snappy_pipeline/config.yaml	2017-02-03T13:12:22.693821
    @@ -33,6 +33,10 @@
         path_ngs_mapping: ../ngs_mapping  # REQUIRED
         scalpel:
         path_target_regions:  # REQUIRED
    +  somatic_variant_annotation:
    +    path_somatic_variant_calling: ../somatic_variant_calling  # REQUIRED
    +    oncotator:
    +      path_corpus: REQUIRED  # REQUIRED
     # Data Sets =======================================================================================
     #
     # Define data sets.  The search paths and patterns are given per data set.
    @@ -50,4 +54,5 @@
     #       - /fast/projects/medgen_genomes/2017-01-09_acheiropodia
     #     type: germline_variants
     #
    -data_sets: {}
    +data_sets                   # REQUIRED
    +: {}
    [...]
    $ tree somatic_variant_annotation
    somatic_variant_annotation
    +-- sge_log/
    `-- config.yaml


Adding Sample Sheets
====================

.. note:: The following does not work yet but should in the future

    **TODO**

For matched cancer studies, the most simple way of creating a sample sheet is starting from the shortcut TSV.
The following creates a sample sheet TSV shortcut.
This is then converted into a JSON bio-med sample sheet.

.. code-block:: shell

    $ cat <<"EOF" | sed $'s/[ \t]\+/\t/g' > .snappy_pipeline/01_data_set.tsv
    [Metadata]
    schema          cancer_matched
    schema_version  v1
    title           Example matched cancer tumor/normal study
    description     The study has two patients, P001 has one tumor sample, P002 has two

    [Data]
    patientName sampleName  isTumor    libraryType folderName
    P001    N1  N   WES P001-N1-DNA1-WES1
    P001    T1  Y   WES P001-T1-DNA1-WES1
    P001    T1  Y   mRNA_seq    P001-T1-RNA1-mRNA_seq1
    P002    N1  N   WES P002-N1-DNA1-WES1
    P002    T1  Y   WES P002-T1-DNA1-WES1
    P002    T1  Y   WES P002-T1-RNA1-mRNA_seq1
    P002    T2  Y   WES P002-T2-DNA1-WES1
    P002    T2  Y   mRNA_seq    P002-T2-RNA1-mRNA_seq1
    EOF
    $ biomedsheets -t matched_cancer \
        --input .snappy_pipeline/01_data_set.tsv \
        --output .snappy_pipeline/01_data_set.json
    $ head .snappy_pipeline/01_data_set.json
    [TODO]

.. note::

    Updating entries in data set TSV files does not work yet and requires a re-starting from scratch.
    As the data set primary keys are part of the file names, changing the PK of sample or library will require cleaning all output files and re-running the whole pipeline.
    Overall, it is better to only use the JSON sheet files and the corresponding tools and helpers.

Now, we have to register the data set in the configuration.
Ensure that the ``data_sets`` entry look as follows.
Replace ``<path-to-demo-dir>`` with the path to the ``demo`` directory of the ``snappy_pipeline`` project.

.. code-block:: yaml

    data_sets:
      first_batch:
        file: 01_first_batch.tsv
        search_patterns:
          # Note that currently only "left" and "right" key known
          - {'left': '*/L???/*_R1.fastq.gz', 'right': '*/L???/*_R2.fastq.gz'}
        search_paths: ['<path-to-demo-dir>/input/01_first_batch']
        type: matched_cancer

.. TODO: describe full configuration and setting format

The full configuration format will be described elsewhere.
It is notable, however, that there also is an optional ``naming_scheme`` property for each batch.
Using this, you can select between naming based on secondary ID and pk (``secondary_id_pk``) and secondary ID alone (``only_secondary_id``).


Executing the Project's Pipeline
================================

After executing the steps from above, our pipeline is ready to use.
Each pipeline step instance will automatically run each predecessor within the pipeline.
Thus, it is enough to execute the pipeline in the ``somatic_variant_annotation`` step.

For running, locally use:

.. code-block:: shell

    $ cd somatic_variant_annotation
    $ cubi-snake -p --step somatic_variant_annotation

For running with SGE on the cluster, use the ``--cubi-pipeline-drmaa`` parameter.

.. code-block:: shell

    $ cd somatic_variant_annotation
    $ cubi-snake -p --step somatic_variant_annotation --cubi-pipeline-drmaa

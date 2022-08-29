===========================
CUBI Pipeline Documentation
===========================

This is the documentatation for the CUBI Pipeline.
This documentation is split into four parts:

Pipeline User Docs
    Documentation for **pipeline users**.
    Start here to learn about the pipeline.
    This section starts at :ref:`quickstart`.

Pipeline Step Docs
    Documentation for the individual pipeline steps.
    This includes a general description, description of the related configuration settings, and a documentation of generated output files and input workflow steps.
    This section starts at :ref:`pipeline_steps_introduction`.

Pipeline Developers Docs
    Documentation for **pipeline developers**.
    After you are proficient in using the pipeline, continue reading here if you want to fix, change, or extend the pipeline.
    This section starts at :ref:`dev_intro`.

API Documentation
    This is the entry point for the API.

Project Info
    House-keeping information about the project, such as instructions for developer setup, author list, changelog etc.
    Start at :ref:`howto_release`).

.. note:: Where to Start?

    Even if you want to modify the pipeline, it's best to read the user documentation first (BIH users start a :ref:`quickstart`, other users refer to :ref:`installation`) as you need to be able to run the pipeline to test your changes and additions.

.. toctree::
    :caption: Pipeline Users Docs
    :hidden:

    quickstart
    installation
    usage
    overview

.. toctree::
    :caption: Pipeline Step Docs
    :hidden:

    step_intro
    step_generic
    step/helper_gcnv_model
    step/hla_typing
    step/igv_session_generation
    step/ngs_data_qc
    step/ngs_mapping
    step/ngs_sanity_checking
    step/roh_calling
    step/somatic_gene_fusion_calling
    step/somatic_neoepitope_prediction
    step/somatic_ngs_sanity_checking
    step/somatic_purity_ploidy_estimate
    step/somatic_targeted_seq_cnv_calling
    step/somatic_variant_annotation
    step/somatic_variant_calling
    step/somatic_variant_checking
    step/somatic_variant_expression
    step/somatic_variant_filtration
    step/somatic_wgs_cnv_calling
    step/somatic_wgs_sv_calling
    step/targeted_seq_cnv_calling
    step/targeted_seq_repeat_analysis
    step/tcell_crg_report
    step/variant_annotation
    step/variant_calling
    step/variant_checking
    step/variant_denovo_filtration
    step/variant_phasing
    step/variant_filtration
    step/wgs_cnv_annotation
    step/wgs_cnv_calling
    step/wgs_mei_calling
    step/wgs_sv_annotation
    step/wgs_sv_calling
    step/wgs_sv_filtration


.. toctree::
    :caption: Pipeline Developers Docs
    :hidden:

    dev_intro
    dev_somatic_variant_calling
    dev_ngs_mapping
    dev_api

.. toctree::
    :caption: Project Info
    :hidden:

    contributing
    howto_release
    authors
    history
    license

.. Generated pages, should not appear

    * :ref:`genindex`
    * :ref:`modindex`
    * :ref:`search`

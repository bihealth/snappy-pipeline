.. _dev_somatic_variant_calling:

==================================
Somatic Variant Calling Dissection
==================================

.. note::

    This chapter was written for the legacy ``somatic_variant_calling`` module.
    The workflow has since been renamed to ``variant_calling`` and its internal
    architecture has changed significantly (single-tool selection, task-based
    config, ``produces``/``consumes`` contracts).

    For the current implementation, see the ``variant_calling`` step source at
    ``snappy_pipeline/workflows/variant_calling/``.

    The general concepts (``BaseStep``, ``BaseStepPart``, input/output file
    patterns) are described in :ref:`dev_intro` and demonstrated in
    :ref:`dev_ngs_mapping`.

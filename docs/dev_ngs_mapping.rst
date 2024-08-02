.. _dev_ngs_mapping:

======================
NGS Mapping Dissection
======================

This chapter gives a dissection of the NGS Mapping step (``ngs_mapping``) for CUBI pipeline developers.
The NGS mapping step is an example for a pipeline step that works on the raw FASTQ NGS read files.
This chapter assumes that you have read :ref:`dev_somatic_variant_calling` before.

The minority of pipeline steps will work directly with the raw read data.
Most steps work on the results of the NGS read mapping step or even further downstream.

.. note::

    Before reading this chapter, you should

    - have knowledge from the user's perspective of CUBI pipeline (start a :ref:`usage`)
    - have read chapter :ref:`dev_intro`
    - have read chapter :ref:`dev_somatic_variant_calling`.

    After reading this chapter, you should

    - know how to work with raw FASTQ file input
    - know how to use the using the :py:class:`LinkInStep <snappy_pipeline.workflows.abstract.LinkInStep>`

**This is still TODO, just look at the code for now ;)**

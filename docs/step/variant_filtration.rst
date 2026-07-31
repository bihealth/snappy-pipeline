.. _step_variant_filtration:

===========================
Germline Variant Filtration
===========================

.. automodule:: snappy_pipeline.workflows.variant_filtration

Overview
========

The ``variant_filtration`` step has been unified to support both germline and somatic workflows under a task-based model. Each configured task executes a single filtration tool. Complex pipelines are built by chaining tasks using the ``depends_on.variant`` parameter.

Supported Tools
===============

The workflow supports the following filtration tools:

- ``bcftools``: Expression-based tagging or filtering.
- ``vembrane``: Python-expression-based tagging (soft filtering) or filtering (hard removal).
- ``regions``: Region/BED-based inclusion or exclusion.
- ``dkfz``: DKFZ bias filter (cancer-specific, requires aligned BAMs).
- ``ebfilter``: EBFilter (cancer-specific, requires aligned BAMs).

Modes: Tagging vs. Filtering
============================

Each tool supports a ``mode`` toggle:

- ``mode: tag`` (Default): Adds a custom filter label to the ``FILTER`` column of VCF records.
- ``mode: filter``: Removes non-matching records completely from the VCF.

Configuration Example
=====================

Here is an example config chaining a quality filter (using vembrane in ``tag`` mode) followed by a frequency and region hard filter:

.. code-block:: yaml

    tasks:
      - step: variant_filtration
        name: quality_tagging
        depends_on:
          variant: variant_annotation_vep
        config:
          tool: vembrane
          vembrane:
            mode: tag
            expressions:
              low_gq: 'any(FORMAT["GQ"][s] < 30 for s in SAMPLES)'
              low_depth: 'any(FORMAT["DP"][s] < 10 for s in SAMPLES)'

      - step: variant_filtration
        name: hard_filtration
        depends_on:
          variant: quality_tagging
        config:
          tool: vembrane
          vembrane:
            mode: filter
            expression: '("low_gq" not in FILTER) and ("low_depth" not in FILTER) and (INFO.get("gnomAD_AF", 0.0) < 0.01)'

Migration Guide: Legacy to Vembrane
===================================

The legacy germline filtration workflow configured complex ``filter_combinations`` from predefined threshold blocks. These can be migrated directly into cleaner and more powerful ``vembrane`` tasks:

1. Quality Thresholds
---------------------

* **Legacy**:
  .. code-block:: yaml

      min_gq: 40
      min_dp_het: 10
      min_dp_hom: 5

* **Vembrane**:

  .. code-block:: python

      # As tag expressions:
      expressions:
        low_gq: 'any(FORMAT["GQ"][s] < 40 for s in SAMPLES)'
        poor_support: 'any((is_het(s) and FORMAT["DP"][s] < 10) or (is_hom(s) and FORMAT["DP"][s] < 5) for s in SAMPLES)'

2. Frequency Thresholds
-----------------------

* **Legacy**:
  .. code-block:: yaml

      af_dominant: 0.001
      ac_dominant: 3

* **Vembrane**:
  .. code-block:: python

      # Check population frequencies and allele counts in INFO:
      expression: 'INFO.get("gnomAD_AF", 0.0) < 0.001 and INFO.get("gnomAD_AC", 0) <= 3'

3. Score Thresholds
-------------------

* **Legacy**:
  .. code-block:: yaml

      require_gerpp_gt2: True
      min_cadd: 15

* **Vembrane**:
  .. code-block:: python

      expression: 'INFO.get("GERP_score", 0.0) > 2.0 and INFO.get("CADD_PHRED", 0.0) >= 15'

4. Coding vs. Non-Coding
------------------------

* **Legacy**:
  .. code-block:: yaml

      require_coding: True

* **Vembrane**:
  .. code-block:: python

      # Inspect VEP/Jannovar ANN field annotations:
      expression: 'any(ann["Consequence"] != "synonymous_variant" for ann in ANN)'



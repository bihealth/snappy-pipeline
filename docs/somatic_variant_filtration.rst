.. _somatic_variant_filtration:

----------------------------
Somatic variation filtration
----------------------------

The complete processing of small variants involves 3 steps:

1. | *Calling*: the tool of choice here is ``mutect2`` (other callers are not maintained anymore), usually with a panel of normals.  
   | The output provides a ``*.full.vcf.gz`` which contains all variants found by ``mutect2``, including those rejected, with the reason for rejection.  
   | Only the variants that pass all the caller's filter are considered for further processing.
2. | *Annotation*: the annotation is usually done with the ENSEMBL Variant Effect Predictor (``vep``).  
   | Here again, the output contains a ``*.full.vcf.gz`` with annotations for all genes & transcripts overlapping with the variant locus.  
   | The reduced output (``*.vcf.gz``) contains only one annotation per variant, for the transcript selected as most representative by ``vep``.
3. | *Filtration*: additional filtration step, fully configurable by the user.  
   | The ``*.full.vcf.gz`` contains all annotated variants, including those rejected, with the reason for rejection.  
   | Only the variants that pass all the filters are considered for furthe processing.

The order of these steps cannot (unfortunately) be altered.

Available filters
=================

There are currently 5 available filters. Each can be used or not, in any order, and can appear multiple times.
The filters are:

- `dkfz <https://github.com/DKFZ-ODCF/DKFZBiasFilter>`_ filters out variants which have severe bias in orientation, either on strand or from read pair.
- `ebfilter <https://github.com/Genomon-Project/EBFilter>`_ uses a Bayesian statistical model to score variants, which can then be rejected on that basis.
- ``bcftools`` can filter variants based on their properties. For example, to filter out variants with a depth less than 50, less than 5 reads supporting the variant, or with a Variant Allele Fraction (VAF) less than 5%, the syntax would be:

  .. code-block:: yaml

      somatic_variant_filtration:
        filter_list:
        - bcftools:
          exclude: "AD[1:0]+AD[1:1]<50 | AD[1:1]<5 | AD[1:1]/(AD[1:0]+AD[1:1])<0.05"

- ``regions`` filters out variants that overlap with a user-defined ``bed`` file. Typically, it can be used to reject variants outside of coding regions.
- ``protected`` is a way to protect variants overlapping with a user-defined ``bed`` file. Typically, such variants are known drivers, or actionable variants, that should appear in the final ``vcf`` file, even if they don't fulfill the (possibly) stringent filtering rules.

The ``bcftools`` is a very powerful and flexible filter, but it requires some care.
The syntax for ``exclude`` & ``include`` configuration options are delegated to the ``bcftools`` program as `filtration expressions <https://samtools.github.io/bcftools/bcftools.html#expressions>`_.
Therefore, the presence of ``vcf`` tags and the order of samples are important and influence the outcome of the filter.
In particular, the ``AD[1:0]`` stands for *number of reads supporting the reference allele in the tumor sample* only if the tumor sample is the second sample in the ``vcf`` file.
This is usually the case when the normal sample is ``donor-N1-DNA1-WES1`` & the tumor sample ``donor-T1-DNA1-WES1``, but it is not always true.

The filtration step will be extended to simplify the common operations.
For example, an additional filter might be added with a simplified syntax, to allow easy filtration based on read depth, support & VAF in the tumor sample.

.. note::

    The filtration method described above is now the preferred way.
    The old method, using ``filter_sets`` is now deprecated, because it is very inflexible and inefficient, as it creates a very large number of files for each sample, with very long filenames.


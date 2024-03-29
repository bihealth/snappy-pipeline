.. _somatic_cnv:

-------------------
Somatic CNV calling
-------------------

Somatic variant calling is implemented differently for exome and whole genome data.

The whole genome data "branch" is currently under review, as GRCh38 support in ``Control-FREEC`` (the main workhorse for WGS CNV calling) is not complete.
CNV calling in WGS data can also be done using ``cnvkit``, but its pipeline implementation is also incomplete.

The following documentation is restricted to the tools currently implemented to process exome data: ``cnvkit``, ``purecn`` & ``sequenza``.

Output & performance
====================

The 3 methods generally broadly agree on the log ratio of coverage between tumor & normal samples. 

However, the segmentation and the number of copies assigned to a segment can be quite different between the algorithms.

Output files
------------

There is no widely used standard to report copy number alterations. 
In absence of a better solution, all CNV tools implemented in somatic pipeline output the segmentation table loosely following the `DNAcopy format <https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf>`_.`
The copy number call may or may not be present, and the chromosome number is replaced by its name.
The segmentation output is in file ``output/<mapper>.<cnv caller>.<lib name>/out/<mapper>.<cnv caller>.<lib name>_dnacopy.seg``.

Genome support
--------------

Both ``purecn`` & ``sequenza`` have better support for ``GRCh37`` than for ``GRCh38``.
In both cases, it stems from the fact that the segmentation packages used by the methods haven not been updated for ``GRCh38``.
Both ``purecn`` & ``sequenza`` provide remedial solutions, however they are not available from ``bioconda``. 
For that reason, the pipeline implementation relies on Docker containers for ``purecn`` and on github forks of R packages for ``sequenza``.

The pipeline currently cannot handle male & female samples differently, even though the sex has a strong influence on the results.
It is therefore useful to avoid calling CNVs on sex chromosomes.
This can be easily achieved for ``sequenza``, using the ``ignore_chroms`` option. 
``purecn`` will handle the sample sex internally, but ``cnvkit`` requires the user to remove the sex chromosomes from the access file.

Segmentation and CNV calls
--------------------------

``cnvkit`` generally produces more segments than ``purecn``, and shorter ones. 
However, with the default settings, ``cnvkit`` is very conservative with the calls.
Most segments are called diploid, even when the coverage log ratio suggests otherwise.
For best results, it is advisable to experiment with different segmentation algorithms and different sensitivity thresholds.

Purity & ploidy
---------------

``purecn`` attempts to model tumor purity and ploidy, and the results are found in ``<mapper>.purecn.<lib name>/out/<mapper>.purecn.<lib name>.csv``.
``sequenza`` can also offer these estimates, but they are currently only stored in ``R`` objects in the ``work`` directory, at ``<mapper>.sequenza.<lib name>/report/<lib name>_sequenza_cp_table.RData``.

Getting results for ``PureCN``
==============================

The pipeline implementation of ``PureCN`` requires a panel of normals (which in turns requires building a panel of normals for ``mutect2``).
Moreover, it also needs ``vcf`` files containing both germline & somatic variants.
For these files, ``mutect2`` needs to be run with different parameters than those used for somatic variant calling.
Therefore, the ``somatic_variant_calling`` step must be run *twice*, with different arguments, in different locations.

For the second ``mutect2`` run executed in the ``somatic_variant_calling_for_purecn`` directory, the configuration would look like:

.. code-block:: yaml

    somatic_mutation_calling:
      tools: [mutect2]
      mutect2:
        extra_arguments: [ "--genotype-germline-sites true", "--genotype-pon-sites true" ] # These arguments must be added
        panel_of_normals: <absolute path to the panel of normal output>
        germline_resource: <path to af-only-gnomad.raw.sites.vcf.gz>
        common_variants: <path to small_exac_common_3.vcf.gz>
        window_length: 300000000       # For exome data, it is sufficient to split the genome by chromosomes
        ignore_chroms: [NC_007605, hs37d5, chrEBV, '*_decoy', 'HLA-*', 'GL000220.*'] # For hs37d5
      ignore_chroms: [NC_007605, hs37d5, chrEBV, '*_decoy', 'HLA-*', 'GL000220.*'] # Must be repeated at the level above mutect2
    
    somatic_targeted_seq_cnv_calling:
      tools: [purecn]
      purecn:
        genome_name: "hg19"             # This must match the names given while building the panel of normals
        enrichment_kit_name: "exome"    # This must match the names given while building the panel of normals
        path_somatic_variants: ../somatic_variant_calling_for_purecn
        somatic_variant_caller: mutect2
        path_panels_of_normals: <absolute path to panel_of_normals/output/<mapper>.purecn/out/<mapper>.purecn.panel_of_normals.rds>
        path_mapping_bias: <absolute path to panel_of_normals/output/<mapper>.purecn/out/<mapper>.purecn.mapping_bias.rds>
        path_intervals: <absolute path to panel_of_normals/output/purecn/out/exome_hg19.list>
        path_container: <absolute path to panel_of_normals/work/containers/out/purecn.simg>

From the ``panel_of_normals`` directory, ``purecn`` requires 3 types of files:

- the ``panel_of_normals`` itself, and the ``mapping_bias`` objects are taken from ``<mapper>.purecn/out``. This is because they might change with different mapping tools.
- the ``intervals`` taken from ``purecn/out``, as the definition of intervals depend only on the genome & the exome kit, but not on the mapping tool.
- the ``container`` taken from ``work/containers/out``, to ensure that the ``PureCN`` version used to compute copy number variants is identical to that used to compute the panel of normals.


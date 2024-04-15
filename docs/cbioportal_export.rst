.. _cbioportal_export:

--------------------
Export to cBioPortal
--------------------

The step ensures that the following operations are carried out:

- For each sample, the somatic variant ``vcf`` file must be converted to the ``maf`` format.
- For each sample, the CNV segments results must be applied to gene locii, for the coverage log ratio, and for the copy number count.
- The individual sample files must then be merged to create single files for the whole cohort. This must be done for somatic variants, coverage log ratio & copy number counts, and for expression levels when present.
- Clinical data must be gathered at the patient & sample level.
- List of cases must be created, for mutations, copy number changes & expression levels.
- Meta files required by cBioPortal must be created.

The ``vcf`` to ``maf`` conversion is done by a bespoke script, which requires a mapping between gene symbols and NCBI (ENTREZ) gene ids.
This mapping is independent of the genome release, and can be obtained directly from `HGNC <https://www.genenames.org/download/archive/>`_. 
There is a relatively recent version on the cluster, in ``/fast/work/groups/cubi/projects/biotools/static_data/annotation/hgnc_complete_set_2023-10-01.tsv``.

Enriching clinical data
=======================

The pipeline can produce several quantities of interest, that can be added to the sample data uploaded to cBioPortal.
Currently, only the tumor mutational burden computed by the ``tumor_mutational_burden`` step can be added to the sample table.
In the near future, the Micro Satellite Instability (MSI) status and the homologous recombination deficiency will be added.
The tumor purity and ploidy might also be considered, although as these results are generally hidden in the output of the somatic copy number steps, the situation is less clear.
Finally, the patient HLA alleles and the sex could be included, either at the patient or sample level.

Functional scores based on expression data will also be included, as well as gene fusions, but structural variants might remain a more distant goal, as they are notoriously diffict to infer from exome data.


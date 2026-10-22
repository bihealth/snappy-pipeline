(cbioportal_export)=

# Export to cBioPortal

The `cbioportal_export` step turns the results of other tasks into a study for cBioPortal.
It does these operations:

- Converts the somatic variant VCF file of each sample to MAF format.
- Applies the CNV segments of each sample to genes, for the coverage log ratio and the copy number count.
- Merges the sample files into one file per data type for the cohort: mutations, log ratios, copy numbers and, if present, expression levels.
- Collects the clinical data of the patients and samples.
- Creates the case lists for mutations, copy number changes and expression levels.
- Creates the meta files that cBioPortal requires.

See {ref}`step_cbioportal_export` for the outputs and all options.

## Inputs

The task reads these tasks through `depends_on`:

`variants`
: Required.
    A task with somatic SNVs and indels annotated by VEP, such as `variant_annotation` with `tool: vep` after `variant_calling`.
    The conversion to MAF cannot read mehari's annotations yet.
    Filtration can come before or after the annotation.

`copy_number`
: A task with somatic copy number segments in DNAcopy format, such as `somatic_targeted_seq_cnv_calling`.
    Needed with `copy_number_alteration.enabled: true`.

`alignments`
: An RNA `ngs_mapping` task with `tool: star`, whose gene counts are the expression levels.
    Needed with `expression.enabled: true`.

`features`
: The gene annotation, an `external_data` task.
    Needed with `copy_number_alteration` or `expression`.

## Configuration

The conversion from VCF to MAF uses a script that needs a mapping from gene symbols to NCBI (Entrez) gene ids.
The mapping does not depend on the genome release.
`hgnc_complete_set.txt` from [HGNC](https://www.genenames.org/download/archive/) provides it, and `path_gene_id_mappings` points to it.

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/genome.fa}

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /path/to/bwa/genome}

  - name: vep_cache
    step: external_data
    config:
      produces: {type: database, tags: [vep]}
      files: {cache: /path/to/vep}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: data_sets
        reference: genome
        index: bwa_index
      tool: bwa
      bwa: {}

  - name: calling
    step: variant_calling
    config:
      depends_on:
        alignments: mapping
        reference: genome
      tool: mutect2
      mutect2:
        contamination: {enabled: false}

  - name: annotation
    step: variant_annotation
    config:
      depends_on:
        variants: calling
        reference: genome
        database: vep_cache
      tool: vep
      vep:
        assembly: GRCh38

  - name: export
    step: cbioportal_export
    config:
      depends_on:
        variants: annotation
      path_gene_id_mappings: /path/to/hgnc_complete_set.txt
      vcf2maf:
        Center: my_center
        ncbi_build: GRCh38
      study:
        type_of_cancer: mixed
        cancer_study_id: mixed_pi_2026
        study_name: My study
        study_name_short: My study
        study_description: Somatic exomes
```

`study.type_of_cancer` is a code from [OncoTree](http://oncotree.mskcc.org/#/home).
`study.cancer_study_id` is usually `<type of cancer>_<PI>_<year>`.

## Clinical data

The clinical data of patients and samples come from the sample sheets.
Results of other tasks, such as the tumor mutational burden, are not yet added to the sample table.
`sample_info` is to be redesigned.

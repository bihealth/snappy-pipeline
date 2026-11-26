(cbioportal_export)=

# Export to cBioPortal

The `cbioportal_export` step turns the results of other tasks into the files of a study that cBioPortal imports.
It does these operations:

- Converts the somatic variant VCF file of each sample to MAF format.
- Converts the CNV segments of each sample into a log ratio of the coverage and a copy number count per gene.
- Merges the sample files into one file per data type for the cohort: mutations, log ratios, copy numbers and, if present, expression levels.
- Collects the clinical data of the patients and samples.
- Creates the case lists for mutations, copy number changes and expression levels.
- Creates the meta files that cBioPortal requires.

See {ref}`step_cbioportal_export` for the outputs and all options.

## Inputs

The task reads other tasks through these `depends_on` keys:

`variants`
: Required.
    A task with somatic SNVs and indels annotated by VEP, such as `variant_annotation` with `tool: vep` after `variant_calling`.
    The conversion to MAF cannot read mehari's annotations yet.
    Filtration can come before or after the annotation.

`copy_number`
: A task with somatic copy number segments in DNAcopy format, such as `somatic_targeted_seq_cnv_calling`.
    Required with `copy_number_alteration.enabled: true`.

`alignments`
: An RNA `ngs_mapping` task with `tool: star`, whose gene counts are the expression levels.
    Required with `expression.enabled: true`.

`features`
: The gene annotation, an `external_data` task.
    Required with `copy_number_alteration` or `expression`.

```{graphviz}
:alt: The tasks that cbioportal_export reads and the depends_on keys that name them

digraph cbioportal {
    rankdir=LR; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;
    label=<<FONT POINT-SIZE="7" COLOR="#898781">Arrow labels: the depends_on key of the reading task, optional ones in parentheses</FONT>>;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    "variant_calling" [label=<<B>variant_calling</B><BR ALIGN="LEFT"/>>];
    "variant_annotation" [label=<<B>variant_annotation</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tool: vep</FONT><BR ALIGN="LEFT"/>>];
    "cnv" [label=<<B>somatic_targeted_seq_cnv_calling</B><BR ALIGN="LEFT"/>>];
    "rna" [label=<<B>ngs_mapping</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">RNA, tool: star</FONT><BR ALIGN="LEFT"/>>];
    "features" [label=<<B>external_data</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">gene annotation</FONT><BR ALIGN="LEFT"/>>];
    "cbioportal_export" [label=<<B>cbioportal_export</B><BR ALIGN="LEFT"/>>];
    "variant_calling" -> "variant_annotation" [color="#1baf7ab3", label="variants"];
    "variant_annotation" -> "cbioportal_export" [color="#1baf7ab3", label="variants"];
    "cnv" -> "cbioportal_export" [color="#1baf7ab3", label="(copy_number)"];
    "rna" -> "cbioportal_export" [color="#eb6834b3", label="(alignments)"];
    "features" -> "cbioportal_export" [color="#2a78d6b3", label="(features)"];
}
```

## Configuration

The copy number and expression tables need a mapping from gene symbols to NCBI (Entrez) gene ids.
The mapping does not depend on the genome release.
`hgnc_complete_set.txt` from [HGNC](https://www.genenames.org/download/archive/) provides it.
`path_gene_id_mappings` points to it.

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/genome.fa}

  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]

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
        reads: reads
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

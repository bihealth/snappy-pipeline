(somatic_ngs)=

# Mapping and friends

Somatic variants in low-purity or degraded samples, and sub-clonal variants, depend on the quality of the mapping.
Adapter trimming and base quality score recalibration (BQSR) can improve it.
snappy handles no molecular barcodes (UMIs) for now.

Each is a task of its own: `adapter_trimming` before the mapping, `alignment_processing` after it.
A chain of tasks combines them, and the variant callers read the BAM files of its last task.

## Adapter trimming

The `bwa`, `bwa_mem2` and `bwa_mem3` sections of `ngs_mapping` have a fast option, `trim_adapters`, to trim adapters on the fly.
Where that is not enough, `adapter_trimming` offers `bbduk` and `fastp`.

Unlike other steps, `adapter_trimming` produces FASTQ files.
A task that reads FASTQ files uses them when its `depends_on.reads` names the `adapter_trimming` task.
These steps read FASTQ files: `ngs_mapping`, `hla_typing`, `ngs_data_qc`, `gene_expression_quantification` and `somatic_gene_fusion_calling`.

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

  - name: trimming
    step: adapter_trimming
    config:
      depends_on:
        reads: data_sets
      tool: bbduk
      bbduk:
        adapter_sequences: [/path/to/adapters.fa]

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: trimming
        reference: genome
        index: bwa_index
      tool: bwa
      bwa: {}
```

The mapping task reads the FASTQ files from `tasks/trimming/work/`.
`adapter_trimming` does not rename the files, so the `search_patterns` of the data set still apply.

## Base quality score recalibration

[BQSR should generally be applied](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).
`alignment_processing` with `tool: gatk4_bqsr` recalibrates the base qualities of the BAM files of a mapping task, or of another `alignment_processing` task.
`depends_on.common_variants` names an `external_data` task with a VCF of common germline variants, for example `small_exac_common_3.vcf` of the GATK resource bundle.

## Recipe: mapping and BQSR

The mapping hands its BAM files on and deletes them once `bqsr` has read them (`delete_intermediates: [bam]`).
Variant callers name `bqsr` in `depends_on.alignments`.

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/genome.fa}

  - name: common_variants
    step: external_data
    config:
      produces: {type: variants, tags: [population, common]}
      files: {vcf: /path/to/small_exac_common_3.vcf.gz}

  - name: bwa_mem2_index
    step: reference_index
    config:
      depends_on: {reference: genome}
      tool: bwa_mem2

  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: data_sets, reference: genome, index: bwa_mem2_index}
      tool: bwa_mem2
      bwa_mem2: {}
      delete_intermediates: [bam]

  - name: bqsr
    step: alignment_processing
    config:
      depends_on: {alignments: mapping, reference: genome, common_variants: common_variants}
      tool: gatk4_bqsr
      gatk4_bqsr: {}
```

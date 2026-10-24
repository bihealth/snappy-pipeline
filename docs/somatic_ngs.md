(somatic_ngs)=

# Mapping and friends

Somatic variants in low-purity or degraded samples, and sub-clonal variants, depend on the quality of the mapping.
Adapter trimming, molecular barcodes (UMIs, MBCs) and base quality score recalibration (BQSR) can improve it.

Each is a task of its own: `adapter_trimming` before the mapping, `alignment_processing` after it.
A chain of tasks combines them, and the variant callers read the BAM files of its last task.

## Adapter trimming

The `bwa`, `bwa_mem2` and `bwa_mem3` sections of `ngs_mapping` have a fast option, `trim_adapters`, to trim adapters on the fly.
Where that is not enough, `adapter_trimming` offers `bbduk` and `fastp`, and `agent` for reads with molecular barcodes.

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

## Barcodes and UMIs

Barcodes are handled in three operations:

1. The barcodes, typically the first bases of a read, are clipped from the read.
   They are added to the read name or to the comment on the identifier line of the read record.
2. Mapping runs as usual.
   Depending on the downstream tool, the barcode sequences and qualities must be added as tags to the BAM file.
3. The aligned reads and the barcode sequences are used for deduplication.
   This can produce consensus reads that are longer than the original reads.

Three tasks do them with Agilent's AGeNT:

1. `adapter_trimming` with `tool: agent` clips the barcodes and moves them into the read names.
2. `ngs_mapping` maps the reads with `extra_arguments: ["-C"]`, so that the barcodes become tags.
3. `alignment_processing` with `tool: agent` marks the duplicates with the barcodes, with the `baits` of each library's kit in `depends_on.kits`.

```{note}
Only [AGeNT](https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/ngs-software/agent-232879) from Agilent is implemented for barcodes.
It is not open source and not in bioconda, so install it separately; the `path` of the AGeNT tasks names its jar files.
Open-source alternatives such as [UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html) may replace it.
The AGeNT documentation comes with its download.
```

## Base quality score recalibration

[BQSR should generally be applied](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).
`alignment_processing` with `tool: gatk4_bqsr` recalibrates the base qualities of the BAM files of a mapping task, or of another `alignment_processing` task.
`depends_on.common_variants` names an `external_data` task with a VCF of common germline variants, for example `small_exac_common_3.vcf` of the GATK resource bundle.

## Recipe: barcodes and BQSR

This chain replaces the `mbcs` tool of `ngs_mapping` in snappy 0.6.
The mapping and the duplicate marking hand their BAM files on and delete them once the next task has read them (`delete_intermediates: [bam]`).
Variant callers name `bqsr` in `depends_on.alignments`.
Without barcodes, the chain is the mapping followed by `bqsr`.

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/genome.fa}

  - name: kits
    step: external_data
    config:
      produces: {type: regions, tags: [kits]}
      kits:
        V5: {targets: /path/to/v5_targets.bed, baits: /path/to/v5_baits.bed}

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

  - name: trimming
    step: adapter_trimming
    config:
      depends_on: {reads: data_sets}
      tool: agent
      agent:
        path: /path/to/agent/trimmer.jar
        lib_prep_type: v2    # halo, hs, xt, v2 or qxt, see the AGeNT documentation
        extra_arguments: ["-polyG 8", "-minFractionRead 50"]

  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: trimming, reference: genome, index: bwa_mem2_index, kits: kits}
      tool: bwa_mem2
      bwa_mem2:
        split_as_secondary: true
        extra_arguments: ["-C"]    # the barcodes become tags
      delete_intermediates: [bam]

  - name: marking
    step: alignment_processing
    config:
      depends_on: {alignments: mapping, reference: genome, kits: kits}
      tool: agent
      agent:
        path: /path/to/agent/creak.jar
        consensus_mode: HYBRID    # SINGLE, DUPLEX or HYBRID
        input_filter_args: ["-mm 13", "-mr 13", "-mq 30"]
      delete_intermediates: [bam]

  - name: bqsr
    step: alignment_processing
    config:
      depends_on: {alignments: marking, reference: genome, common_variants: common_variants}
      tool: gatk4_bqsr
      gatk4_bqsr: {}
```

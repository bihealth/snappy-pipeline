(somatic_ngs)=

# Mapping and friends

Somatic variants in low-purity or degraded samples, and sub-clonal variants, depend on the quality of the mapping.
Adapter trimming, molecular barcodes (UMIs, MBCs) and base quality score recalibration (BQSR) can improve it.

Adapter trimming is its own step, `adapter_trimming`.
Barcodes and BQSR are part of the `ngs_mapping` step, in the `mbcs` tool.

## Adapter trimming

The `bwa` and `bwa_mem2` sections of `ngs_mapping` have a fast option, `trim_adapters`, to trim adapters on the fly.
Where that is not enough, `adapter_trimming` offers two tools, `bbduk` and `fastp`.

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

The `mbcs` tool does all three operations in one job per library, followed by optional BQSR.
`mbcs.mapping_tool` selects the mapper, `bwa` or `bwa_mem2`.
The pipeline is only tested with these short-read mappers.
The mapper section (`bwa` or `bwa_mem2`) configures the mapping.

The task takes the index of its `mapping_tool`: `depends_on.index` names a task tagged with that tool, such as `bwa_mem2`.

```{note}
Only [AGeNT](https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/ngs-software/agent-232879) from Agilent is implemented for barcodes.
It is not open source and not in bioconda, so install it separately.
Open-source alternatives such as [UMI-tools](https://umi-tools.readthedocs.io/en/latest/index.html) may replace it.
```

An `mbcs` task with AGeNT looks like this.
The `agent` paths are inputs of the job.

```yaml
- step: ngs_mapping
  name: mapping
  config:
    depends_on:
      reads: data_sets
      reference: genome
      index: bwa_mem2_index
      kits: kits    # external_data task with the targets and baits of each kit
      common_variants: common_variants    # external_data task with the VCF of common germline variants, for BQSR
    tool: mbcs
    mbcs:
      mapping_tool: bwa_mem2
      barcode_tool: agent    # the only barcode tool
      use_barcodes: true
      recalibrate: true
    bwa_mem2:
      trim_adapters: false
      split_as_secondary: true
      extra_args: ["-C"]     # -C when AGeNT processed barcodes are present
    agent:
      prepare:
        path: /path/to/agent/trimmer
        lib_prep_type: v2    # halo, hs, xt, v2 or qxt, see the AGeNT documentation
        extra_args:
          - "-polyG 8"
          - "-minFractionRead 50"
      mark_duplicates:
        path: /path/to/agent/creak
        consensus_mode: HYBRID    # SINGLE, DUPLEX or HYBRID
        input_filter_args: ["-mm 13", "-mr 13", "-mq 30"]
        consensus_filter_args: []
        extra_args: []
```

The AGeNT documentation comes with its download.
The AGeNT baits file is the `baits` file of the library's kit in the `kits` task.
Each library names its kit in the sample sheet column `libraryKit`.

## Base quality score recalibration

[BQSR should generally be applied](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).
It is implemented only in the `mbcs` tool, with `recalibrate: true` (the default).
Without barcodes, set `use_barcodes: false` (the default).
`depends_on.common_variants` names an `external_data` task with a VCF of common germline variants, for example `small_exac_common_3.vcf` of the GATK resource bundle.

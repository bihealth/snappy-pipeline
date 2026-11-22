(somatic_ngs)=

# Mapping and friends

Somatic variants in low-purity or degraded samples depend on the quality of the mapping.
So do sub-clonal variants.
Adapter trimming and base quality score recalibration (BQSR) can improve the quality of the mapping.
snappy handles no molecular barcodes (UMIs) for now.

Both are tasks of their own: `adapter_trimming` before the mapping, `alignment_processing` after it.
A chain of tasks combines them.
The variant callers read the BAM files of the last task of the chain.

```{graphviz}
:alt: The chain of tasks from the FASTQ files of the data sets to the variant callers

digraph chain {
    rankdir=LR; nodesep=0.12; ranksep=0.3;
    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;
    label=<<FONT POINT-SIZE="7" COLOR="#898781">Arrow labels: the depends_on key of the reading task</FONT>>;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    "data_sets" [style="rounded,dashed", label=<<B>data_sets</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">FASTQ files of the sample sheet</FONT><BR ALIGN="LEFT"/>>];
    "adapter_trimming" [label=<<B>adapter_trimming</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">optional, tool: bbduk or fastp</FONT><BR ALIGN="LEFT"/>>];
    "ngs_mapping" [label=<<B>ngs_mapping</B><BR ALIGN="LEFT"/>>];
    "alignment_processing" [label=<<B>alignment_processing</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">optional, tool: gatk4_bqsr</FONT><BR ALIGN="LEFT"/>>];
    "callers" [label=<<B>variant callers</B><BR ALIGN="LEFT"/>>];
    "data_sets" -> "adapter_trimming" [color="#2a78d6b3", label="reads"];
    "adapter_trimming" -> "ngs_mapping" [color="#2a78d6b3", label="reads"];
    "ngs_mapping" -> "alignment_processing" [color="#eb6834b3", label="alignments"];
    "alignment_processing" -> "callers" [color="#eb6834b3", label="alignments"];
}
```

## Adapter trimming

`trim_adapters` in the `bwa`, `bwa_mem2` and `bwa_mem3` sections of `ngs_mapping` trims adapters on the fly.
It is a fast option.
Where that is not enough, `adapter_trimming` offers `bbduk` and `fastp`.

Unlike other steps, `adapter_trimming` produces FASTQ files.
A task that reads FASTQ files uses them when its `depends_on.reads` names the `adapter_trimming` task.
These steps read FASTQ files:

- `ngs_mapping`
- `hla_typing`
- `ngs_data_qc`
- `gene_expression_quantification`
- `somatic_gene_fusion_calling`

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
`alignment_processing` with `tool: gatk4_bqsr` recalibrates the base qualities of BAM files.
It reads them from a mapping task or from another `alignment_processing` task.
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

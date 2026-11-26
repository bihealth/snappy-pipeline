(somatic_ngs)=

# Mapping and friends

Somatic variants in low-purity or degraded samples depend on the quality of the mapping.
So do sub-clonal variants.
Adapter trimming and base quality score recalibration (BQSR) can improve the quality of the mapping.
Libraries with unique molecular identifiers (UMIs) take an `umi_mapping` task instead of `ngs_mapping` (see below).

Adapter trimming and BQSR are tasks of their own: `adapter_trimming` runs before the mapping, `alignment_processing` after it.
Each task names the task before it in `depends_on`, so together they form a chain.
The variant callers read the BAM files of the last task of the chain.

```{graphviz}
:alt: The chain of tasks from the FASTQ files of the libraries to the variant callers

digraph chain {
    rankdir=LR; nodesep=0.12; ranksep=0.3;
    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;
    label=<<FONT POINT-SIZE="7" COLOR="#898781">Arrow labels: the depends_on key of the reading task</FONT>>;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    "reads" [style="rounded,dashed", label=<<B>reads</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">external_data: the FASTQ files of the libraries</FONT><BR ALIGN="LEFT"/>>];
    "adapter_trimming" [label=<<B>adapter_trimming</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">optional, tool: bbduk or fastp</FONT><BR ALIGN="LEFT"/>>];
    "ngs_mapping" [label=<<B>ngs_mapping</B><BR ALIGN="LEFT"/>>];
    "umi_mapping" [label=<<B>umi_mapping</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">UMI libraries, tool: fgumi</FONT><BR ALIGN="LEFT"/>>];
    "alignment_processing" [label=<<B>alignment_processing</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">optional, tool: gatk4_bqsr</FONT><BR ALIGN="LEFT"/>>];
    "callers" [label=<<B>variant callers</B><BR ALIGN="LEFT"/>>];
    "reads" -> "adapter_trimming" [color="#2a78d6b3", label="reads"];
    "adapter_trimming" -> "ngs_mapping" [color="#2a78d6b3", label="reads"];
    "ngs_mapping" -> "alignment_processing" [color="#eb6834b3", label="alignments"];
    "reads" -> "umi_mapping" [color="#2a78d6b3", label="reads"];
    "umi_mapping" -> "callers" [color="#eb6834b3", label="alignments"];
    "alignment_processing" -> "callers" [color="#eb6834b3", label="alignments"];
}
```

## Adapter trimming

`trim_adapters` in the `bwa`, `bwa_mem2` and `bwa_mem3` sections of `ngs_mapping` trims adapters on the fly, in the mapping job.
This is fast and writes no trimmed FASTQ files.
Where it is not enough, an `adapter_trimming` task trims with `bbduk` or `fastp`.

Unlike other steps, `adapter_trimming` produces FASTQ files.
A task reads the trimmed files when its `depends_on.reads` names the `adapter_trimming` task.
These steps read FASTQ files:

- `ngs_mapping`
- `umi_mapping`
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

  - name: reads
    step: external_data
    config:
      produces: {type: raw, tags: [fastq]}
      library_files: [reads.tsv]
      scan:
        search_paths: [raw]
        search_patterns:
          - left: '(?P<readgroup>.+)\.R1\.fastq\.gz'
            right: '(?P<readgroup>.+)\.R2\.fastq\.gz'

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /path/to/bwa/genome}

  - name: trimming
    step: adapter_trimming
    config:
      depends_on:
        reads: reads
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

The mapping task reads the trimmed FASTQ files from `tasks/trimming/work/`.
`adapter_trimming` does not rename the files, so the read groups of the files table still apply.

## Base quality score recalibration

[BQSR should generally be applied](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR).
`alignment_processing` with `tool: gatk4_bqsr` recalibrates the base qualities of BAM files.
It reads the BAM files of the task in `depends_on.alignments`: a mapping task or another `alignment_processing` task.
`depends_on.common_variants` names an `external_data` task with a VCF of common germline variants, for example `small_exac_common_3.vcf` of the GATK resource bundle.
BQSR does not count these known variant sites as sequencing errors.

## UMI consensus mapping

Libraries with unique molecular identifiers (UMIs) are mapped by `umi_mapping` with `tool: fgumi` instead of `ngs_mapping` ({ref}`step_umi_mapping`).
The task extracts the UMIs, aligns the reads with bwa-mem3, groups them by UMI and position, calls one consensus read per molecule, and realigns and filters the consensus reads.
Its result files are BAM files like those of a mapping task, so the variant callers and the other consumers read them through `depends_on.alignments` as usual.
`consensus: duplex` serves libraries with UMIs on both reads, whose two strands form one duplex consensus read; `consensus: simplex` serves libraries with one UMI per molecule.
The read structures of the FASTQ files come from the sample sheet's custom field `readStructure`, for example `8M+T 8M+T` for 8 bp UMIs at the start of both reads.
A kit with a fixed set of UMIs gets a `correct` section and an `external_data` task with the UMIs in `depends_on.umis`.

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
      scan:
        search_paths: [raw]
        search_patterns:
          - left: '(?P<readgroup>.+)\.R1\.fastq\.gz'
            right: '(?P<readgroup>.+)\.R2\.fastq\.gz'

  - name: umis
    step: external_data
    config:
      produces: {type: tabular, tags: [umis]}
      files: {umis: /path/to/kit_umis.txt}

  - name: bwa_mem3_index
    step: reference_index
    config:
      depends_on: {reference: genome}
      tool: bwa_mem3

  - name: mapping
    step: umi_mapping
    config:
      depends_on: {reads: reads, reference: genome, index: bwa_mem3_index, umis: umis}
      tool: fgumi
      consensus: duplex
      fgumi:
        correct: {max_mismatches: 1, min_distance: 2}
        duplex: {}
```

## Recipe: mapping and BQSR

With `delete_intermediates: [bam]`, the mapping task deletes its BAM files once `bqsr` has read them.
Variant callers name `bqsr` in `depends_on.alignments`.

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
      scan:
        search_paths: [raw]
        search_patterns:
          - left: '(?P<readgroup>.+)\.R1\.fastq\.gz'
            right: '(?P<readgroup>.+)\.R2\.fastq\.gz'

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
      depends_on: {reads: reads, reference: genome, index: bwa_mem2_index}
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

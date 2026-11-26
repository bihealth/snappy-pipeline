(panel_of_normals)=

# Panel of normals

A panel of normals uses normal samples to find technical biases in the data.
Later steps correct these biases in the tumor samples.
Collect, process and sequence all normal samples of a panel in the same way.
Use the panel only for tumor samples that were collected, processed and sequenced like the panel.

A panel needs enough normal samples.
With too few normals, features of single normals count as technical bias, and later steps correct the tumor samples for them.

## Panel tasks

The `panel_of_normals` step creates panels for somatic variant calling (`mutect2`) and for somatic copy number calling of exome data (`cnvkit`, `purecn`).
Each panel is its own task, with one `tool`.
A calling task names its panel task in `depends_on.panel_of_normals`.
See {ref}`step_panel_of_normals` for the outputs and all options.

```{graphviz}
:alt: The panel of normals tasks of this page, what they read and which tasks read them

digraph panels {
    rankdir=LR; newrank=true; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    "mapping" [label=<<B>ngs_mapping</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">BAM files of the normals</FONT><BR ALIGN="LEFT"/>>];
    "pon_access" [label=<<B>pon_access</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: access</FONT><BR ALIGN="LEFT"/>>];
    "pon_mutect2" [label=<<B>pon_mutect2</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: mutect2</FONT><BR ALIGN="LEFT"/>>];
    "pon_cnvkit" [label=<<B>pon_cnvkit</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: cnvkit</FONT><BR ALIGN="LEFT"/>>];
    "pon_purecn" [label=<<B>pon_purecn</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: purecn</FONT><BR ALIGN="LEFT"/>>];
    "variant_calling" [label=<<B>variant_calling</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tool: mutect2</FONT><BR ALIGN="LEFT"/>>];
    "targeted" [label=<<B>somatic_targeted_seq_cnv_calling</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tool: cnvkit or purecn</FONT><BR ALIGN="LEFT"/>>];
    "wgs" [label=<<B>somatic_wgs_cnv_calling</B><BR ALIGN="LEFT"/>>];
    { rank=same; "mapping"; "pon_access"; }
    { rank=same; "variant_calling"; "targeted"; "wgs"; }
    "mapping" -> "pon_mutect2" [color="#eb6834b3"];
    "mapping" -> "pon_purecn" [color="#eb6834b3"];
    "mapping" -> "pon_cnvkit" [color="#eb6834b3"];
    "pon_access" -> "pon_cnvkit" [color="#2a78d6b3"];
    "pon_mutect2" -> "variant_calling" [color="#008300b3"];
    "pon_mutect2" -> "pon_purecn" [color="#008300b3", label="GenomicsDB"];
    "pon_purecn" -> "targeted" [color="#008300b3"];
    "pon_cnvkit" -> "targeted" [color="#008300b3", label="exome"];
    "pon_cnvkit" -> "wgs" [color="#008300b3", label="whole genome"];
}
```

## Somatic variant calling

The `mutect2` panel needs a germline resource in `depends_on.germline_resource`: a VCF with the population allele frequencies of known germline variants.
The [GATK best practice bundle](https://console.cloud.google.com/storage/browser/gatk-best-practices) provides it.
Look for VCF files with `af-only-gnomad` in their name.

GATK [recommends](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON) at least 40 normals.
Smaller panels are common.
For projects with too few normals, the same bundle has a panel of normals from GATK.

By default, the panel is built from all normal DNA libraries of the sample sheets, those of background data sets included.
Two options change what goes into it:

`path_normals_list`
: A file with one library name per line.
    The panel is built from these libraries instead, even if the sample sheet does not flag them as normal.

`ignore_contigs`
: Regular expressions, matched with `re.search` against the contig names.
    The panel leaves out the matching contigs and so covers only the sequences of interest.

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

  - name: gnomad
    step: external_data
    config:
      produces: {type: variants, tags: [population]}
      files: {vcf: /path/to/af-only-gnomad.raw.sites.vcf.gz}

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /path/to/bwa/genome}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: reads
        reference: genome
        index: bwa_index
      tool: bwa
      bwa: {}

  - name: pon_mutect2
    step: panel_of_normals
    config:
      depends_on:
        alignments: mapping
        reference: genome
        germline_resource: gnomad
      tool: mutect2
      ignore_contigs: ['^NC_007605$', '^hs37d5$', '^chrEBV$', '_decoy$', '^HLA-', '^GL000220']  # for hs37d5
      mutect2:
        path_normals_list: ""   # empty: all normal libraries
```

The panel is `tasks/pon_mutect2/work/panel_of_normals/out/panel_of_normals.vcf.gz`.
A `variant_calling` task with `tool: mutect2` reads it when its `depends_on.panel_of_normals` names the panel task.

## Copy numbers of exome data

The panels for `cnvkit` and `purecn` need extra files, such as the regions of the exome kit.
Their calling tasks take these files from the panel task and do not name them again.

### cnvkit

Without `depends_on.kits`, `cnvkit` runs in whole genome mode.
For exomes, `depends_on.kits` names an `external_data` task with `kits`.
The targets are then the baits of the normals' kit.
All normals must be targeted libraries of the same kit.

A `somatic_targeted_seq_cnv_calling` task reads an exome panel in `depends_on.panel_of_normals`.
A `somatic_wgs_cnv_calling` task reads a whole genome panel in the same key.

An access file lists the regions that are suitable for coverage analysis, without repeats, hard-to-map regions or extreme GC content.
CNVkit provides one for `GRCh37`, but it is incomplete for `GRCh38`.
A task with `tool: access` creates it from the reference and from BED files of regions to exclude (`access.exclude`), such as mappability regions.
It writes the access file to `tasks/<access task>/work/panel_of_normals/out/panel_of_normals.accessible.bed`.

The `cnvkit` task names the access task in `depends_on.access`.
An `external_data` task with `produces: {type: raw, tags: [access]}` and `files: {access: <BED file>}` can take its place.
The access file is optional.
Without it, exome mode does not restrict the antitarget regions.
Whole genome mode creates an access file with default settings instead.

```yaml
- name: pon_access
  step: panel_of_normals
  config:
    depends_on:
      reference: genome
    tool: access
    access:
      exclude: [/path/to/mappability.bed]

- name: pon_cnvkit
  step: panel_of_normals
  config:
    depends_on:
      alignments: mapping
      reference: genome
      access: pon_access
      kits: kits   # omit for whole genome data
    tool: cnvkit
    cnvkit: {}
```

The panel task writes the reports `panel_of_normals.metrics.tsv` and `panel_of_normals.sex.tsv` to `panel_of_normals/report/`.
The CNVkit authors recommend checking them before the panel is used.
The pipeline cannot yet create separate panels for male and female samples.

### PureCN

`purecn` reads the GenomicsDB of a `mutect2` panel, so it needs a `mutect2` panel task.
It also needs `depends_on.kits`, because the baits of the normals' kit are its bait regions.
Create the `mutect2` task first, and name it in `depends_on.panel_of_normals`:

```yaml
- name: pon_purecn
  step: panel_of_normals
  config:
    depends_on:
      alignments: mapping
      reference: genome
      panel_of_normals: pon_mutect2
      kits: kits
    tool: purecn
    purecn:
      genome_name: hg19             # hg19 or hg38
      enrichment_kit_name: exome    # used to name output files
```

PureCN runs in a container image, because its R packages are too complex for a conda environment.
Snakemake pulls the image `purecn.container` and runs the PureCN jobs inside it, see {ref}`usage`.
The default image is pinned by its digest.
`purecn.container` also takes the path of a local `.sif` file.
See {ref}`somatic_cnv` for how a calling task uses the panel.

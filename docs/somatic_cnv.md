(somatic_cnv)=

# Somatic CNV calling

Somatic copy number calling is implemented separately for targeted sequencing (exomes and panels) and for whole genomes.
This page covers the targeted tools of `somatic_targeted_seq_cnv_calling`: `cnvkit`, `purecn` and `sequenza`.
The whole genome step is `somatic_wgs_cnv_calling`, see its step page.

## How the tools differ

All three tools generally agree on the log ratio of the coverage between tumor and normal.
The segmentation and the number of copies assigned to a segment can differ considerably between them.

### Output

There is no widely used standard to report copy number alterations.
All tools write the segmentation table loosely following the [DNAcopy format](https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf).
The copy number call may be missing, and chromosomes are named instead of numbered.
The segmentation of a tumor library is `tasks/<task>/output/<library>/out/<library>.dnacopy.seg`.

### Genome support

`purecn` and `sequenza` support `GRCh37` better than `GRCh38`.
The segmentation packages they use have not been updated for `GRCh38`.
Both offer remedies that are not in bioconda.
So `purecn` runs in a container, and the environment of `sequenza` installs the R package from its repository at a fixed commit when it is created.

### Sex chromosomes

The pipeline cannot yet treat male and female samples differently, although the sex strongly influences the results.
Avoid calling CNVs on the sex chromosomes:

- `sequenza` has the `ignore_contigs` option.
- `purecn` handles the sex of a sample itself.
- `cnvkit` needs the sex chromosomes removed from the access file.

### Segmentation and calls

`cnvkit` usually finds more and shorter segments than `purecn`.
With default settings, `cnvkit` calls conservatively, and most segments are called diploid even if the log ratio suggests otherwise.
Try different segmentation algorithms and sensitivity thresholds for best results.

### Purity and ploidy

`purecn` models tumor purity and ploidy.
They are in `tasks/<task>/output/<library>/out/<library>.csv`.

`sequenza` estimates them too.
They are only stored as R objects in `tasks/<task>/work/<library>/report/<library>_sequenza_cp_table.RData`.

## Configuration

`cnvkit` and `purecn` need a panel of normals, see {ref}`panel_of_normals`.

### cnvkit

A `cnvkit` task names the `cnvkit` panel task in `depends_on.panel_of_normals`.
It reads the target and antitarget regions from the panel task, so the config has no paths to them:

```yaml
- name: cnv_cnvkit
  step: somatic_targeted_seq_cnv_calling
  config:
    depends_on:
      alignments: mapping
      reference: genome
      panel_of_normals: pon_cnvkit
    tool: cnvkit
    cnvkit: {}
```

### PureCN

`purecn` needs three things:

1. A `purecn` panel of normals, which needs a `mutect2` panel of normals.
2. VCF files with germline and somatic variants.
   Mutect2 must call these with other options than the usual somatic calling.
   So the `variant_calling` step runs a second time, as its own task.
3. The container image of the panel task (`purecn.container`), so that the panel and the calls use the same PureCN.
   The `cnv_purecn` task runs its PureCN jobs in the image of the task in `depends_on.panel_of_normals`.

The second `mutect2` task needs `--genotype-germline-sites true` and `--genotype-pon-sites true`, and the `mutect2` panel of normals.

```{graphviz}
:alt: The tasks of the PureCN example and the depends_on keys that connect them

digraph purecn {
    rankdir=LR; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;
    label=<<FONT POINT-SIZE="7" COLOR="#898781">Arrow labels: the depends_on key of the reading task; unlabeled arrows: alignments; external_data tasks not shown</FONT>>;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    "mapping" [label=<<B>mapping</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">ngs_mapping</FONT><BR ALIGN="LEFT"/>>];
    "pon_mutect2" [label=<<B>pon_mutect2</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: mutect2</FONT><BR ALIGN="LEFT"/>>];
    "pon_purecn" [label=<<B>pon_purecn</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">panel_of_normals, tool: purecn</FONT><BR ALIGN="LEFT"/>>];
    "calling_for_purecn" [label=<<B>calling_for_purecn</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">variant_calling, tool: mutect2</FONT><BR ALIGN="LEFT"/>>];
    "cnv_purecn" [label=<<B>cnv_purecn</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">somatic_targeted_seq_cnv_calling, tool: purecn</FONT><BR ALIGN="LEFT"/>>];
    "mapping" -> "pon_mutect2" [color="#eb6834b3"];
    "mapping" -> "pon_purecn" [color="#eb6834b3"];
    "mapping" -> "calling_for_purecn" [color="#eb6834b3"];
    "mapping" -> "cnv_purecn" [color="#eb6834b3"];
    "pon_mutect2" -> "pon_purecn" [color="#008300b3", label="panel_of_normals"];
    "pon_mutect2" -> "calling_for_purecn" [color="#008300b3", label="panel_of_normals"];
    "pon_purecn" -> "cnv_purecn" [color="#008300b3", label="panel_of_normals"];
    "calling_for_purecn" -> "cnv_purecn" [color="#1baf7ab3", label="variants"];
}
```

```yaml
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/genome.fa}

  - name: gnomad
    step: external_data
    config:
      produces: {type: variants, tags: [population]}
      files: {vcf: /path/to/af-only-gnomad.raw.sites.vcf.gz}

  - name: common_variants
    step: external_data
    config:
      produces: {type: variants, tags: [common]}
      files: {vcf: /path/to/small_exac_common_3.vcf.gz}

  - name: kits
    step: external_data
    config:
      produces: {type: regions, tags: [kits]}
      kits:
        MedExome: {targets: /path/to/exome_targets.bed, baits: /path/to/exome_baits.bed}

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /path/to/bwa/genome}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on:
        reads: data_sets
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
      mutect2: {}

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
        genome_name: hg19
        enrichment_kit_name: exome

  - name: calling_for_purecn
    step: variant_calling
    config:
      depends_on:
        alignments: mapping
        reference: genome
        panel_of_normals: pon_mutect2
        germline_resource: gnomad
        common_variants: common_variants
      tool: mutect2
      mutect2:
        extra_arguments:
          - "--genotype-germline-sites true"
          - "--genotype-pon-sites true"
        contamination:
          enabled: true

  - name: cnv_purecn
    step: somatic_targeted_seq_cnv_calling
    config:
      depends_on:
        alignments: mapping
        variants: calling_for_purecn
        panel_of_normals: pon_purecn
      tool: purecn
      purecn: {}
```

The `cnv_purecn` task takes these from `pon_purecn`, so that the panel and the calls agree:

- the normal database, the mapping bias and the intervals
- the container image, in which it runs PureCN
- the genome and the seed of the panel task's `purecn` section

The intervals depend only on the genome and the exome kit.
The normal database and the mapping bias can change with the mapper.

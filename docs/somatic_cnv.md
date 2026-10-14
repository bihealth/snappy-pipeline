(somatic_cnv)=

# Somatic CNV calling

Somatic copy number calling is implemented separately for targeted sequencing (exomes and panels) and for whole genomes.
This page covers the targeted tools of `somatic_targeted_seq_cnv_calling`: `cnvkit`, `purecn` and `sequenza`.
The whole genome step is `somatic_wgs_cnv_calling`, see its step page.

All three tools generally agree on the log ratio of the coverage between tumor and normal.
The segmentation and the number of copies assigned to a segment can differ considerably between them.

## Output

There is no widely used standard to report copy number alterations.
All tools write the segmentation table loosely following the [DNAcopy format](https://bioconductor.org/packages/devel/bioc/manuals/DNAcopy/man/DNAcopy.pdf).
The copy number call may be missing, and chromosomes are named instead of numbered.
The segmentation of a tumor library is `tasks/<task>/output/<library>/out/<library>.dnacopy.seg`.

## Genome support

`purecn` and `sequenza` support `GRCh37` better than `GRCh38`.
The segmentation packages they use have not been updated for `GRCh38`.
Both offer remedies that are not in bioconda.
So `purecn` runs in a container, and `sequenza` uses GitHub forks of R packages.

The pipeline cannot yet treat male and female samples differently, although the sex strongly influences the results.
Avoid calling CNVs on the sex chromosomes:

- `sequenza` has the `ignore_chroms` option.
- `purecn` handles the sex of a sample itself.
- `cnvkit` needs the sex chromosomes removed from the access file.

## Segmentation and calls

`cnvkit` usually finds more and shorter segments than `purecn`.
With default settings, `cnvkit` calls conservatively, and most segments are called diploid even if the log ratio suggests otherwise.
Try different segmentation algorithms and sensitivity thresholds for best results.

## Purity and ploidy

`purecn` models tumor purity and ploidy.
They are in `tasks/<task>/output/<library>/out/<library>.csv`.

`sequenza` estimates them too.
They are only stored as R objects in `tasks/<task>/work/<library>/report/<library>_sequenza_cp_table.RData`.

## cnvkit

`cnvkit` and `purecn` need a panel of normals, see {ref}`panel_of_normals`.
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

## PureCN

`purecn` needs three things:

1. A `purecn` panel of normals, which needs a `mutect2` panel of normals.
2. VCF files with germline and somatic variants.
   Mutect2 must call these with other options than the usual somatic calling, so the `variant_calling` step runs a second time as its own task.
3. The container of the panel task, so that the PureCN version of the panel and the calls are the same.
   The `cnv_purecn` task reads it from `depends_on.panel_of_normals`.

The second `mutect2` task needs `--genotype-germline-sites true` and `--genotype-pon-sites true`, and the `mutect2` panel of normals.

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
      produces: {type: raw, tags: [kits]}
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
      tool: mutect2
      mutect2:
        germline_resource: /path/to/af-only-gnomad.raw.sites.vcf.gz

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
      tool: mutect2
      mutect2:
        extra_arguments:
          - "--genotype-germline-sites true"
          - "--genotype-pon-sites true"
        germline_resource: /path/to/af-only-gnomad.raw.sites.vcf.gz
        contamination:
          enabled: true
          common_variants: /path/to/small_exac_common_3.vcf.gz

  - name: cnv_purecn
    step: somatic_targeted_seq_cnv_calling
    config:
      depends_on:
        alignments: mapping
        reference: genome
        variants: calling_for_purecn
        panel_of_normals: pon_purecn
      tool: purecn
      purecn:
        genome_name: hg19
        enrichment_kit_name: exome
        somatic_variant_caller: mutect2
```

The `cnv_purecn` task reads the normal database, the mapping bias, the intervals and the container from `pon_purecn`.
The intervals depend only on the genome and the exome kit.
The normal database and the mapping bias can change with the mapper.

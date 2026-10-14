(panel_of_normals)=

# Panel of normals

A panel of normals uses normal samples to find technical biases in the data.
Later steps correct these biases in the tumor samples.
Collect, process and sequence all normal samples of a panel in the same way.
Use the panel only for tumor samples that were collected, processed and sequenced like the panel.

A panel needs enough normal samples.
With too few, the particulars of single normals count as bias and are corrected in the tumor samples.

The `panel_of_normals` step creates panels for somatic variant calling (`mutect2`) and for somatic copy number calling of exome data (`cnvkit`, `purecn`).
Each panel is its own task, with one `tool`.
A copy number task names the panel task in `depends_on.panel_of_normals`.

See {ref}`step_panel_of_normals` for the outputs and all options.

## Somatic variant calling

The panel for `mutect2` needs a germline resource, a VCF with the allele frequencies of known germline variants.
The [GATK best practice bundle](https://console.cloud.google.com/storage/browser/gatk-best-practices) provides it.
Look for VCF files with `af-only-gnomad` in their name.

GATK [recommends](https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON) at least 40 normals.
Smaller panels are common.
With too few samples, GATK provides a panel of normals in the same bundle.

By default, all normal DNA libraries of the sample sheets build the panel.
`path_normals_list` selects other libraries: a file with one library name per line.
These are used even if the sample sheet does not flag them as normal.

`ignore_chroms` restricts the panel to the sequences of interest.
It lists regular expressions, matched with `re.search` against the contig names.

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
      ignore_chroms: ['^NC_007605$', '^hs37d5$', '^chrEBV$', '_decoy$', '^HLA-', '^GL000220']  # for hs37d5
      mutect2:
        germline_resource: /path/to/af-only-gnomad.raw.sites.vcf.gz
        path_normals_list: ""   # empty: all normal libraries
```

The panel is `tasks/pon_mutect2/work/panel_of_normals/out/panel_of_normals.vcf.gz`.
A `mutect2` variant calling task takes it from the panel task in `depends_on.panel_of_normals`.

## Copy numbers of exome data

The panels for `cnvkit` and `purecn` need extra files, and their calling tasks read these files from the panel task.

### cnvkit

Without `depends_on.kits`, `cnvkit` runs in whole genome mode.
For exomes, `depends_on.kits` names an `external_data` task with `kits`, and the baits of the kit of the normals are the targets.
All normals must be targeted libraries of the same kit.

An access file lists the regions that are suitable for coverage analysis, without repeats, hard-to-map regions or extreme GC content.
CNVkit provides one for `GRCh37` but not for `GRCh38`.
A task with `tool: access` creates it from the reference and from BED files of regions to exclude, such as mappability regions.
Its file is `tasks/<access task>/work/panel_of_normals/out/panel_of_normals.accessible.bed`.
The `cnvkit` task names the access task in `depends_on.access`.
An `external_data` task with `produces: {type: raw, tags: [access]}` and `files: {access: <BED file>}` serves as well.

The access file is optional.
Without it, exome mode does not restrict the antitarget regions, and whole genome mode creates an access file with default settings.

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

The panel task writes `metrics.tsv` and `sex.tsv` reports in `panel_of_normals/report/`.
The CNVkit authors recommend checking them before the panel is used.
The pipeline cannot yet create separate panels for male and female samples.

### PureCN

`purecn` needs the `mutect2` panel, because it uses the GenomicsDB of the `mutect2` panel.
It also needs `depends_on.kits`: the baits of the kit of the normals are its bait regions.
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

PureCN runs in a container, because its R packages are too complex for a conda environment.
The panel task installs the container in `work/containers/out/container.simg`.
See {ref}`somatic_cnv` for how a calling task uses the panel.

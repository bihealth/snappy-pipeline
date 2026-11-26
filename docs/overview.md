(overview)=

# Overview

This chapter explains the concepts of the pipeline system.
It assumes that you know a workflow system such as GNU Make, Snakemake or bpipe, and its benefit over shell scripts.

Automation in bioinformatics pays off as long as the same steps repeat in every project.
Beyond that point, the work becomes project-specific, and shared code helps less.
Report generation is an example, because most of its code cannot be reused.
This pipeline system covers the steps before that point.

(motivation)=

## Motivation

The pipeline system aims at these properties:

Re-use
: Common analysis steps are reused.
    Most of them are calls of standard tools with some glue and conversion code.

Configurability
: No path is hard-coded.
    Paths and the parameters that may need tuning come from the configuration.

Sensible defaults
: The configuration has defaults, ideally tuned from the data.

Documentation
: Good documentation pays off because the steps are widely reused.

Logging software versions
: The versions of the pipeline and the tools are known, so that an analysis can be repeated.
    The versions also explain slight differences in results.

Stable result paths
: The paths in `output/` do not change or disappear between minor versions.

Robustness
: Failures are prevented where possible, for example by checking the configuration before the first job starts.
    Failed jobs are restarted with more memory and time.

Restartability
: A restart or new input data does not repeat finished work.
    A changed intermediate file updates the files that depend on it, as in GNU Make.

Ease of use
: Local and cluster execution are simple.
    Resource requests have defaults.

(definitions)=

## Definitions

pipeline
: Code for a set of bioinformatics analyses, run automatically.

project
: A directory with a `config.yaml`.
    A project is an **instance** of the pipeline: its configuration plugs steps together.

(pipeline) step
: The code (Snakefiles, wrappers, Python classes) for one self-contained stage of an analysis, such as read mapping or variant calling.
    `snappy task list --available` lists the steps.

task
: A configured use of a step in a project, with a unique name.
    A project can have several tasks of the same step, for example two variant callers.

tool
: The program that a task runs.
    A task has exactly one `tool`, such as `bwa` for `ngs_mapping`.
    To use two tools, define two tasks.

(pipeline_projects)=

## Components of a project

A project uses these pieces:

- The `snappy_pipeline` Python package, installed once.
  It provides the `snappy` command and the code of all steps.
- `config.yaml` in the project directory.
  Its `tasks` list defines the tasks, and `data_sets` points to the sample sheets and the FASTQ files.
  The `project` section is required and holds the name, a description and the `external_id`.
  `external_id` identifies the project in a data management system that tracks where the raw data lives, for example as the SODAR project UUID.
- The sample sheets, which describe the libraries and the experimental design.
- The raw data, for example a directory with the FASTQ files of each library.
- Reference data, such as the genome, the gene annotation and the mapper indices.
  `external_data` tasks name existing files, and `reference_download` and `reference_index` tasks create them.
- `tasks/<task_name>/`, which `snappy run` creates for each task, with `work/` and `output/`.

```{graphviz}
digraph {
    rankdir=LR; bgcolor="transparent"; nodesep=0.2; ranksep=0.5;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b"];
    edge [color="#898781", arrowsize=0.6, fontname="Helvetica", fontsize=8, fontcolor="#52514e"];
    config [label="config.yaml\nproject, data_sets, tasks"];
    sheets [label="sample sheets\nlibraries, pedigrees, tumor-normal pairs"];
    raw [label="raw data\nFASTQ files below the search paths"];
    refs [label="reference data\nexternal_data, reference_download, reference_index tasks"];
    run [label="snappy run", style="rounded,filled", fillcolor="#e1e0d9"];
    work [label="tasks/<task>/work/\nall files of a task; other tasks read here"];
    output [label="tasks/<task>/output/\nlinks to the result files, .manifest.tsv"];
    config -> run; sheets -> run; raw -> run; refs -> run;
    run -> work; work -> output [label="results job"];
}
```

`snappy run` validates the configuration once when it starts.
It checks the task configurations, the `depends_on` links and dependency cycles, and whether each upstream task provides what its consumer needs.

## Tasks and dependencies

A task names the tasks it reads from in `config.depends_on`.
The keys are named after the data: `reads`, `alignments`, `variants`, `reference`, `index`, `features`, `dbsnp`, `panel_of_normals` and so on.
Each step documents the keys it accepts, see {ref}`pipeline_steps_introduction`.
A `depends_on` key that the task's tool does not read is an error, and so is a key that only an off setting reads.
A config therefore never names a task it does not use.
A `depends_on` key has no default task name.

A task declares what it `produces`: a data type with tags, such as alignments of DNA.
A consumer checks this when the project loads, so a mapping task with `tool: star` cannot feed a DNA variant caller.

In this small somatic project, mapping reads the FASTQ files of the `data_sets`, and calling reads the alignments:

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
```

The pipelines in `.tests/test-workflow/pipelines/` are complete examples with longer chains (mapping, calling, annotation, filtration).

### External data

Existing data from outside the project comes from `external_data` tasks.
Examples are reads, BAM or VCF files of another pipeline, and project-wide files such as the reference genome.
Such a task has no rules, and consumers read the files where they are.

With `files: {key: path}`, a task provides project-wide files.
With `library_files`, it provides the files of each library from files tables: one row per library, the columns `library` and the keys under which consumers read the files (`vcf`, `bam`, ...), or one row per read group for reads.
`snappy scan-files` writes a table from a `scan` section, as for reads.

```yaml
tasks:
  - name: dragen_calls
    step: external_data
    config:
      produces: {type: variants, tags: [germline, snv, indel]}
      library_files: [dragen_calls.tsv]   # the columns library, vcf, vcf_tbi
      scan:
        search_paths: [/path/to/dragen]
        search_patterns:
          - { vcf: '.+\.hard-filtered\.vcf\.gz', vcf_tbi: '.+\.hard-filtered\.vcf\.gz\.tbi' }
```

### Dependencies and task-local files

A task gets its input files in one of two ways.
A dependency comes from another task, named in `depends_on`.
A task-local file is a path in the task's tool section.

A dependency is data that the pipeline knows.
It has a type and tags, and the consumer checks both, and the keys it reads, when the project loads.
Reference data is a dependency too: the genome, the gene annotation, enrichment kits, exons and other regions, and the indexes of the mappers.
`external_data` tasks provide existing files, and `reference_index` tasks build indexes.
A task that reads an index from a `reference_index` task must name the genome of that index in `depends_on.reference`.
Otherwise the project does not load.
A custom genome, such as mouse with extra contigs or human with bacterial genomes, is one FASTA file that both tasks name.
A genome FASTA needs its index and sequence dictionary next to it, named as samtools and GATK name them (`genome.fa.fai`, `genome.dict`).
Every job that reads the genome has both as inputs, so a missing one fails when Snakemake builds the DAG.
`reference_download` builds both for the FASTA it downloads.
In the same way, a VCF of dbSNP or population variants needs its index next to it, as in the GATK resource bundle.
The index is `.tbi` for a `.vcf.gz` and `.idx` for a plain `.vcf`.
The BED file of mappable regions that gCNV reads needs its index too (`.tbi` for a `.bed.gz`, `.idx` from `gatk IndexFeatureFile` for a plain `.bed`).
Several tasks can read one `external_data` or `reference_index` task.
So each reference file is named once, and a switch to another genome build only changes these tasks.

A file that only one tool reads can still be a dependency.
A mapper index has a format that only its mapper reads, but it must match the genome that the other tasks use.
Its type is `index`, and its tags name the tool and the molecule, as for `bwa_index` in the example above.

A task-local file belongs to the setup of one tool and does not depend on the reference or other data of the project.
Examples are a jar or a Python file with helper functions for vembrane.
Snappy passes the file to the rule as an input.
Snakemake checks that it exists when a job needs it and reruns jobs when it changes, but nothing checks what the file is.

To decide, ask whether the file must match the reference or other data of the project.
If so, it is a dependency, even if only one tool reads it.
Otherwise it is a task-local file.

```{graphviz}
digraph {
    rankdir=LR; bgcolor="transparent"; nodesep=0.3; ranksep=0.5;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b"];
    edge [color="#898781", arrowsize=0.6, fontname="Helvetica", fontsize=8, fontcolor="#52514e"];
    question [label="Must the file match the reference\nor other data of the project?", shape=diamond, style="filled", fillcolor="#e1e0d9"];
    dependency [label="dependency\nan external_data task, named in depends_on\ntype, tags and keys are checked"];
    local [label="task-local file\na path in the tool section\npassed to the rule as an input"];
    question -> dependency [label="yes"];
    question -> local [label="no"];
}
```

BED files of regions have the type `regions` and the key `bed`.
The regions filter of `variant_filtration` accepts any of them in `depends_on.regions`.
Steps that need a certain kind of regions require a tag, for example `exons` for the strandedness of STAR.

This project tags the variants in blacklisted regions.
It then keeps the variants for which a function from a local Python file returns true:

```yaml
tasks:
  - name: calls
    step: external_data
    config:
      produces: {type: variants, tags: [germline, snv, indel]}
      library_files: [calls.tsv]   # the columns library, vcf, vcf_tbi

  - name: blacklist
    step: external_data
    config:
      produces: {type: regions, tags: [blacklist]}
      files: {bed: /path/to/blacklist.bed}

  - name: blacklist_filter
    step: variant_filtration
    config:
      depends_on:
        variants: calls
        regions: blacklist  # a dependency
      tool: regions
      regions: {action: exclude, mode: tag}

  - name: lab_filter
    step: variant_filtration
    config:
      depends_on:
        variants: blacklist_filter
      tool: vembrane
      vembrane:
        mode: filter
        expression: 'is_reportable(QUAL, INFO["DP"])'
        context_files: [filters/lab_rules.py]  # a task-local file
```

(config_paths)=

### Paths

Relative paths in the configuration are relative to the directory of `config.yaml`.
This holds for every path: sample sheets, search paths, `external_data` files and files in tool sections, also when Snakemake runs in another directory.
Absolute paths are not changed.

The top-level `pathvars` section names roots for paths, used as `<name>` in any path:

```yaml
pathvars:
  refs: /path/to/refs/GRCh38
data_sets: {}
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: <refs>/genome.fa}
```

The section works like Snakemake's `pathvars`, with its default names: `<resources>` stands for `resources` unless the section sets it.
An undefined name fails when the project loads.

snappy does not check whether the files exist when the project loads, so a dry-run works on a machine without them.
Snakemake checks the inputs of the jobs it runs.
`snappy check-files` lists every file of the configuration that does not exist, the files of the files tables included.
This includes the `.fai` and `.dict` of a genome and the index of a VCF or of the mappability BED file.

## How the files of the libraries are found

Every task that reads FASTQ files names its source in `depends_on.reads`:

- An `external_data` task that produces `raw` data, with the FASTQ files of the libraries.
- An `adapter_trimming` task, for its trimmed files.

An `external_data` task lists the files of each library in files tables, `library_files`.
A files table is a TSV file with one row per read group of a library, usually a lane: the columns `library`, `readgroup`, `left` and, for paired-end data, `right`.
The `readgroup` becomes the read group of the alignments.
An empty cell is no value, and relative paths are relative to the table's directory.
Read groups are processed in the order of the path of their left file.
Every read group needs a right file when the table has a `right` column, unless `mixed_se_pe: true` on the task accepts read groups without one, so that single-end and paired-end data can be mixed.
Steps read the files where they are, so the project holds no copies or links of FASTQ files.

`snappy scan-files` writes the tables.
The task's `scan` section says where to look: the files of a library are in a folder named like its `folder_name` in the sample sheet, else like the library, at any depth below a `search_paths` entry.
`search_patterns` are regular expressions for the path below that folder, and each must match the whole path.
The `(?P<readgroup>...)` group pairs the `left` and `right` file of a read group.
The command searches once and lists what it found; the table is the record of the run, so keep it in git.
Run it again when the sheets gain libraries.
A library without rows stops `snappy run` and `snappy check-files` with an error that names the libraries and the command; `snappy run --scan` runs the scan first.
A row for a library that no sheet has is an error as well.
Two libraries that share a folder are an error of the scan.

```yaml
data_sets:
  my_batch:
    file: samplesheet.tsv
    type: matched_cancer

tasks:
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

  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /refs/GRCh38.fa}

  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /refs/bwa/GRCh38}

  - name: mapping
    step: ngs_mapping
    config:
      depends_on: {reads: reads, reference: genome, index: bwa_index}
      tool: bwa
      bwa: {}
```

The same tables serve the other entry points: an `external_data` task with `produces: {type: alignments, tags: [dna]}` lists one row per library with the columns `library`, `bam` and `bai`, and a task with `type: variants` the columns `library`, `vcf` and `vcf_tbi`.
Such a project needs no FASTQ files.

## Sample sheets

Sample sheets are TSV files in the [BioMed Sample Sheet](https://github.com/bihealth/biomedsheets) format.
The `type` of the data set selects how the sheet is read: `matched_cancer`, `germline_variants` or `generic`.

A `matched_cancer` sheet has one normal and one tumor sample per patient:

```tsv
[Metadata]
schema          cancer_matched
schema_version  v1

[Data]
patientName  sampleName  isTumor  libraryType  folderName
P001         N1          N        WES          P001-N1-DNA1-WES1
P001         T1          Y        WES          P001-T1-DNA1-WES1
P002         N1          N        WES          P002-N1-DNA1-WES1
P002         T1          Y        WES          P002-T1-DNA1-WES1
```

In a `germline_variants` or `generic` sheet, a row may leave `libraryType` empty.
Such a row describes a person without sequencing data, such as a parent who completes a pedigree.
A sheet in which no row has a `libraryType` has no libraries, and snappy stops with an error.

### Library selection

Each step takes only libraries it can process, such as RNA libraries for expression quantification or tumor DNA libraries for somatic calling.
The somatic steps take only libraries of cancer sample sheets.

Without `library_selection`, a task takes the step's default:

- in general, the tumor DNA libraries of cancer sample sheets and the DNA libraries of germline sample sheets.
- another set where the step's page names one.
  For example, mapping takes all libraries, the RNA steps take all RNA libraries, and a panel of normals takes all normal DNA libraries.

`library_selection` replaces that default with a pandas `DataFrame.query()` expression, for example `data_set == 'batch2'`, `donor_name != 'P003'` or `batch_no == 'run-2024-A'`.
The custom fields of the sheet are columns too, in snake case, with the values of the sheet: a `batchNo` field declared as `string` keeps batch ids, and `umi_mapping` reads the read structures of a UMI library from `readStructure`.
A task never takes libraries that its step cannot process.

### Group by

`group_by` sets what one output covers: `library` writes one output per library, `related` one per group of related libraries, where the step supports it.
Related libraries are the members of a pedigree on germline sample sheets and the libraries of one donor on cancer sample sheets (the library table's column `related_group`).

The group is named after the pedigree or donor, and so are the results of `group_by: related` tasks:

- the sheet's family id (the data set's `pedigree_field`, else the `familyId` custom field),
- else the index donor, with `FAM_` on germline sheets, `PAT_` on cancer sheets and no prefix on generic sheets.

Without `group_by`, a task takes its step's default, else what its upstream tasks declare, else `library`.
Germline variant calling defaults to `related` and calls each pedigree jointly, so the tasks after it also write one result per pedigree.
Steps that read one VCF per pedigree (VarFish export, peddy, IGV sessions, pedigree phasing) check when the project loads that their upstream task has `group_by: related`.

### Background data sets

A data set with `is_background: true` holds libraries that are no task's results, such as a reference cohort.
Mapping and panels of normals process these libraries, and the SV and gCNV callers add them as cohort members or model samples.
`background_selection` narrows them per task, for example to `data_set == 'cohort_1000g'`.

Several data sets may have libraries of the same name.
A task names its files by library, so it can take only one library of a name.
A task and its upstream tasks must take that library from the same data set.
If a task would take both, it fails and names the libraries.
Select one data set with `library_selection`, for example `data_set == 'batch1'`, or use one task per data set.

## Files of a task

Each task has its own directory `tasks/<task_name>/`:

```text
my_project/
+-- config.yaml
+-- tasks/
    +-- mapping/
    |   +-- work/
    |   +-- output/
    +-- calling/
        +-- work/
        +-- output/
```

All files of a task are in `work/`, and other tasks read them there.
`output/` only holds symlinks to the result files, and `output/.manifest.tsv`.
See {ref}`usage` and {ref}`generic_pipeline_step`.

## Resource requests

Each step declares the threads, memory and runtime of its jobs.
`restart-times` in the SLURM profile sets how often a failed job is retried.
Memory and runtime grow with each retry by `resources.retry_factor`:

```yaml
resources:
  retry_factor: 1.5   # default
```

Partition limits are left to SLURM and its executor plugin.
Every job that runs a wrapper writes a benchmark (`tasks/<task>/benchmarks/...tsv`), except jobs of rules cached between workflows.
The benchmark holds the job's wall time (`s`), peak memory (`max_rss`), I/O and CPU time.
It records what each job used, so the requests can be checked against it.

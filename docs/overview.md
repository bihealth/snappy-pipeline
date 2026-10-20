(overview)=

# Overview

This chapter describes the concepts of the pipeline system.
It is for readers who know workflow systems such as GNU Make, Snakemake or bpipe, and who see their benefit over shell scripts.

Automation in bioinformatics pays off while the same steps repeat in every project.
Beyond that point, the work becomes project-specific and shared code helps less.
Report generation is an example, where most code cannot be reused.
This pipeline system covers the steps before that point.

(motivation)=

## Motivation

The pipeline system aims at these properties:

Re-use
: Common analysis steps are reused.
    Mostly these are calls of standard tools with some glue and conversion code.

Configurability
: No path is hard-coded.
    Paths and the parameters that may need tuning come from the configuration.

Sensible defaults
: The configuration has defaults, ideally tuned from the data.

Documentation
: Wide reuse repays good documentation.

Logging software versions
: The versions of the pipeline and the tools are known, so that an analysis can be repeated.
    They also explain slight differences in results.

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
: Code for a set of bioinformatics tasks, run automatically.

project
: A directory with a `config.yaml`.
    A project is an **instance** of the pipeline: the configuration plugs steps together.

(pipeline) step
: The code (Snakefiles, wrappers, Python classes) for one encapsulated set of operations, such as read mapping or variant calling.
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

A project is a directory.
The pieces are:

- The `snappy_pipeline` Python package, installed once.
  It provides the `snappy` command and the code of all steps.
- `config.yaml` in the project directory.
  Its `tasks` list defines the tasks, and `data_sets` points to the sample sheets and the FASTQ files.
- The sample sheets, which describe the libraries and the experimental design.
- The raw data, for example a directory with the FASTQ files of each library.
- Reference data, such as the genome, the gene annotation and the mapper indices.
  `external_data` tasks name existing files, and `reference_download` and `reference_index` tasks create them.
- `tasks/<task_name>/`, which `snappy run` creates for each task, with `work/` and `output/`.

The configuration is validated once when `snappy run` starts.
It checks the task configurations, the `depends_on` links, dependency cycles and that each upstream task provides what its consumer needs.

## Tasks and dependencies

A task names the tasks it reads from in `config.depends_on`.
The keys are named after the data: `reads`, `alignments`, `variants`, `reference`, `index`, `features`, `dbsnp`, `panel_of_normals` and so on.
Each step documents the keys it accepts, see {ref}`pipeline_steps_introduction`.
There are no defaults for task names.

A task declares what it `produces`: a data type with tags, such as alignments of DNA.
A consumer checks this when the project loads, so a mapping task with `tool: star` cannot feed a DNA variant caller.

This is a small somatic project.
Mapping reads the FASTQ files of the `data_sets`, and calling reads the alignments:

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
With `search_paths` and `search_patterns`, it provides per-library files from a folder named like the library.
The pattern keys (`vcf`, `bam`, `left`, ...) are the keys under which consumers read the files.

```yaml
tasks:
  - name: dragen_calls
    step: external_data
    config:
      produces: {type: variants, tags: [germline, snv, indel]}
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
A genome FASTA needs its index and sequence dictionary next to it, named as samtools and GATK name them (`genome.fa.fai`, `genome.dict`).
Every job that reads the genome has both as inputs, so a missing one fails when Snakemake builds the DAG; `reference_download` builds them for the FASTA it downloads.
Several tasks can read one such task.
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
      search_paths: [/path/to/calls]
      search_patterns:
        - { vcf: '.+\.vcf\.gz', vcf_tbi: '.+\.vcf\.gz\.tbi' }

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

Relative paths in the configuration are relative to the project directory, the directory of `config.yaml`.
This holds for every path: sample sheets, search paths, `external_data` files and files in tool sections, also when Snakemake runs in another directory.
Absolute paths are not changed.

The top-level `pathvars` section names roots for paths, used as `<name>` in any path:

```yaml
pathvars:
  refs: /path/to/refs/GRCh38
tasks:
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: <refs>/genome.fa}
data_sets: {}
```

It has the defaults and semantics of Snakemake's `pathvars`, so `<resources>` stands for `resources` unless the section sets it.
An undefined name fails when the project loads.

snappy does not check whether the files exist when the project loads, so a dry-run works on a machine without them.
Snakemake checks the inputs of the jobs it runs.
`snappy check-files` lists every file of the configuration that does not exist, including the `.fai` and `.dict` of a genome.

## How FASTQ files are found

Every task that reads FASTQ files names its source in `depends_on.reads`:

- `data_sets`, for the reads of the libraries in the sample sheets.
- An `external_data` task that produces `raw` data, for reads from elsewhere.
- An `adapter_trimming` task, for its trimmed files.

A data set and an `external_data` task find reads the same way, with the same fields: `search_paths`, `search_patterns` and `mixed_se_pe`.
Each library has a folder, which can be at any depth below a search path.
For a data set, the `folderName` column of the sample sheet names the folder, and only the search paths of the library's own data set are searched.
For an `external_data` task, the folder is named like the library.
The files are found once per `snappy run`.

`search_patterns` are regular expressions for the path below the library folder, and each must match the whole path.
The `(?P<readgroup>...)` group pairs the `left` and `right` file of a read group, usually a lane.
Read groups are processed in the order of the path of their left file.
Every left file needs a right file, unless `mixed_se_pe: true` accepts read groups without one, so that single-end and paired-end data can be mixed.
A data set without `search_patterns` finds `*.R1.fastq.gz` and `*.R2.fastq.gz` (`(?P<readgroup>.+)\.R1\.fastq\.gz` and `...R2...`); `search_paths` have no default.
Steps read the files directly, so there are no copies or links of FASTQ files in the project.

```yaml
data_sets:
  my_batch:
    file: samplesheet.tsv
    type: matched_cancer
    search_paths:
      - raw
    search_patterns:
      - left: '(?P<readgroup>.+)\.R1\.fastq\.gz'
        right: '(?P<readgroup>.+)\.R2\.fastq\.gz'
```

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

Each step takes only libraries it can process, such as RNA libraries for expression quantification or tumor DNA libraries for somatic calling.
The somatic steps take only libraries of cancer sample sheets.
Without `library_selection`, a task takes the step's default: the tumor DNA libraries of cancer sample sheets and the DNA libraries of germline sample sheets, unless the step's page names another one.
For example, mapping takes all libraries, the RNA steps all RNA libraries, and a panel of normals all normal DNA libraries.
`library_selection` replaces that default with a pandas `DataFrame.query()` expression, for example `data_set == 'batch2'` or `donor_name != 'P003'`, but a task never takes libraries that its step cannot process.
`group_by` sets the granularity: `library` writes one output per library, `cohort` one per cohort (pedigree), where the step supports it.
Without it, a task takes its step's default, else the granularity of its upstream tasks, else `library`.
Germline variant calling defaults to `cohort` and calls each pedigree jointly, so the tasks after it also write per cohort.
Steps that read one VCF per pedigree (VarFish export, peddy, IGV sessions, pedigree phasing) check when the project loads that their upstream task writes per cohort.

Several data sets may have libraries of the same name.
A task names its files by library, so it can only take one library of a name, and a task and its upstream tasks must take it from the same data set.
If a task would take both, it fails and names the libraries; select one data set with `library_selection`, for example `data_set == 'batch1'`, or use one task per data set.

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
`output/` only holds symlinks to the published files, and `output/.manifest.tsv`.
See {ref}`usage` and {ref}`generic_pipeline_step`.

## Resource requests

Each step declares the threads, memory and runtime of its jobs.
Memory and runtime grow with each retry of a failed job, by `resources.retry_factor` (`restart-times` in the SLURM profile):

```yaml
resources:
  retry_factor: 1.5   # default
```

Partition limits are left to SLURM and its executor plugin.
Every job that runs a wrapper writes a benchmark (`tasks/<task>/benchmarks/...tsv`) with its wall time (`s`), peak memory (`max_rss`), I/O and CPU time.
These measurements calibrate the requests.

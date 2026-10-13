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

### Relative paths

Relative paths in the configuration are looked up in the project directory, then in its parent directory.
Absolute paths are not changed.

## How FASTQ files are found

Every task that reads FASTQ files names its source in `depends_on.reads`:

- `data_sets`, the folders below the `search_paths` of the data sets.
- An `external_data` task that produces `raw` data.
- An `adapter_trimming` task, for its trimmed files.

The files are found once per `snappy run`.
For `data_sets`, the library folder is named like the `folderName` of the sample sheet, or like the library if the sheet has no `folderName`.
It can be at any depth below a `search_paths` entry.

`search_patterns` are regular expressions for the path below the library folder, and each must match the whole path.
The `(?P<readgroup>...)` group pairs the `left` and `right` file of a read group, usually a lane.
Read groups are processed in the order of the path of their left file.
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

A task processes all libraries of the data sets by default.
`library_selection` restricts a task to the libraries that match a pandas `DataFrame.query()` expression.
`group_by: cohort` makes a task write one output per cohort (pedigree) instead of one per library, where the step supports it.

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

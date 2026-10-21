(quickstart)=

# Quickstart

This chapter sets up a first project: a matched tumor/normal exome, mapped with BWA and called with Mutect2.
It assumes that you have a reference genome with a BWA index, a sample sheet and the FASTQ files.

```{note}
This describes the setup as a pipeline user.
For the developer setup, see {ref}`installation`.
```

## Install Snappy Pipeline

Install [pixi](https://pixi.sh) (see <https://pixi.sh/latest/#installation>).
Then clone the repository and let pixi create the environment:

```shell
git clone git@github.com:bihealth/snappy-pipeline.git
cd snappy-pipeline
pixi install
```

The `snappy` command is then available as `pixi run snappy <subcommand> ...`.
`pixi install --frozen` installs the exact versions in `pixi.lock`.
See {ref}`installation` for details.

## Create the project

`snappy init` creates the project directory with `config.yaml`, `pipeline_job.sh` and `README.md`:

```shell
snappy init my_project
cd my_project
```

Add one task per step with `snappy task add`.
Each argument is `[task_name=]step_type`.
`snappy task list --available` lists the step types.

```shell
snappy task add mapping=ngs_mapping calling=variant_calling
```

The new tasks appear in `config.yaml` with every option commented out, and the required ones marked `REQUIRED`.
`--config-mode full` writes all defaults, `--config-mode minimal` strips them.
`snappy task add` and `snappy init --task` keep the previous `config.yaml` as `config.yaml.bak`.

## Configure the tasks

A project has one `config.yaml` with two sections.
`tasks` lists the tasks, each with a `step`, a unique `name` and a `config`.
`data_sets` names the sample sheets and where their FASTQ files are.

A task names the tasks it reads from in `config.depends_on`.
The keys are named after the data: `reads`, `alignments`, `variants`, `reference`, `index` and so on.
Each value is the name of a task, except `reads: data_sets`, which reads the FASTQ files of the `data_sets`.
Every task sets one `tool` and configures it in the section of the same name.

Reference data comes from `external_data` tasks.
They have no rules and declare what the files are in `produces`.

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
      files: {index: /path/to/bwa/genome}  # prefix of the .amb, .ann, ... files

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

Relative paths are relative to the project directory, the directory of `config.yaml`; `pathvars` can name roots for paths, see {ref}`config_paths`.
The genome FASTA needs its index and sequence dictionary next to it (`genome.fa.fai`, `genome.dict`).
The index can also come from a `reference_index` task, and the reference from a `reference_download` task; a `reference_index` task builds the index from the genome that the mapping task names in `depends_on.reference`.

## Add the sample sheet and the reads

The sample sheet is a TSV file in the [BioMed Sample Sheet](https://github.com/bihealth/biomedsheets) format.
This one is a `matched_cancer` sheet:

```tsv
[Metadata]
schema          cancer_matched
schema_version  v1

[Data]
patientName  sampleName  isTumor  libraryType  folderName
P001         N1          N        WES          P001-N1-DNA1-WES1
P001         T1          Y        WES          P001-T1-DNA1-WES1
```

The reads of a library are in a folder named like its `folderName`, at any depth below a `search_paths` entry of its data set.
`search_patterns` are regular expressions for the path below that folder.
The `(?P<readgroup>...)` group pairs the `left` and `right` file of a read group, usually a lane.

```text
my_project/
+-- config.yaml
+-- samplesheet.tsv
+-- raw/
    +-- P001-N1-DNA1-WES1/
    |   +-- lane1.R1.fastq.gz
    |   +-- lane1.R2.fastq.gz
    +-- P001-T1-DNA1-WES1/
        +-- lane1.R1.fastq.gz
        +-- lane1.R2.fastq.gz
```

`snappy pull-sheet` fetches a sample sheet from SODAR.

## Run

`snappy check-files` lists the files that `config.yaml` names but that do not exist.
Then start with a dry-run.
It validates `config.yaml` and lists the jobs without running them.
Everything after `--` goes to Snakemake.

```shell
snappy run -- --dry-run --cores 1
```

Then run locally, or submit through SLURM:

```shell
snappy run -- --cores 8
snappy run --slurm
```

`pipeline_job.sh` wraps `snappy run --slurm` for `sbatch`.
Results of a task are in `tasks/<task_name>/work/`.
`tasks/<task_name>/output/` links the published files, see {ref}`usage`.

(quickstart)=

# Quickstart

This chapter sets up a first project: a matched tumor/normal exome, mapped with BWA and called with Mutect2.
You need a reference genome with a BWA index, a sample sheet and the FASTQ files.

```{note}
This chapter covers the setup for pipeline users.
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

`snappy init` creates the project directory with `config.yaml`, `README.md` and an empty sample sheet `samplesheet.tsv`.
`config.yaml` names this sheet in `sample_sheets` and a `reads` task for its FASTQ files below `raw/`:

```shell
snappy init my_project
cd my_project
```

On a terminal, `snappy init` asks for what no option gives:

- the project name, its external id in a data management system such as the SODAR project UUID, and a description,
- the sample sheet type (`generic`, `germline` or `cancer`),
- the tasks,
- the git settings.

`--no-interactive` takes the defaults instead.

The project is a git repository whose first commit holds `.gitignore`, `config.yaml` and `README.md`.
Everything else is ignored: results, logs, caches, data, the sample sheet and the files tables, which may hold personal data.
`--track-sheets` adds the sample sheet and the files table `reads.tsv`, so that the record of a run covers the files of the libraries.
Use it only if they hold no personal data or the remote is approved for such data.
`snappy scan-files` points out a table that the repository misses.
`--remote <url>` adds the remote `origin` without pushing, and `--no-git` skips the repository.
Inside an existing git repository, `snappy init` stages the files there instead of creating a new repository.

Add a task for each step with `snappy task add`.
Each argument is `[task_name=]step_type`.
`snappy task list --available` lists the step types.

```shell
snappy task add mapping=ngs_mapping calling=variant_calling
```

The new tasks appear in `config.yaml` with their required fields, marked `REQUIRED`.
`--config-mode commented` writes every field with its description, the optional ones commented out, and `--config-mode full` writes every field active.
Keep the project in git for the history of its configuration.

## Configure the tasks

Besides `project`, `config.yaml` has two main sections.
`sample_sheets` names the sample sheets, by default `samplesheet.tsv` next to `config.yaml`.
`tasks` lists the tasks, each with a `step`, a unique `name` and a `config`.

A task names the tasks it reads from in `config.depends_on`.
The keys are named after the data: `reads`, `alignments`, `variants`, `reference`, `index` and so on.
Each value is the name of a task.
Every task sets one `tool` and configures it in the section of the same name.

Reference data and the files of the libraries come from `external_data` tasks.
They have no rules and declare in `produces` what their files are.
The `reads` task below lists the FASTQ files of each library in `reads.tsv`, which `snappy scan-files` writes.

```yaml
project:
  name: my_project

sample_sheets: [samplesheet.tsv]

tasks:
  - name: reads
    step: external_data
    config:
      produces: {type: reads}
      library_files: [reads.tsv]

  - name: genome
    step: external_data
    config:
      produces: {type: sequences, tags: [reference, dna]}
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

The tasks and what they read from each other:

```{graphviz}
digraph {
    rankdir=LR; bgcolor="transparent"; nodesep=0.2; ranksep=0.5;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b"];
    edge [color="#898781", arrowsize=0.6, fontname="Helvetica", fontsize=8, fontcolor="#52514e"];
    reads [label="reads\nexternal_data, FASTQ files"];
    genome [label="genome\nexternal_data"];
    bwa_index [label="bwa_index\nexternal_data"];
    mapping [label="mapping\nngs_mapping, bwa"];
    calling [label="calling\nvariant_calling, mutect2"];
    reads -> mapping [label="reads"];
    genome -> mapping [label="reference"];
    bwa_index -> mapping [label="index"];
    mapping -> calling [label="alignments"];
    genome -> calling [label="reference"];
}
```

Relative paths are relative to the directory of `config.yaml`.
`pathvars` can name roots for paths, see {ref}`config_paths`.
The genome FASTA needs its index and sequence dictionary next to it (`genome.fa.fai`, `genome.dict`).
A `reference_index` task can build the index instead, and a `reference_download` task can provide the genome.
The mapping task must then name the genome of that index in `depends_on.reference`.

## Add the sample sheet and the reads

The sample sheet is a plain TSV file with one row per library; `snappy init` writes its header.
The core columns are `donor_name`, `library_name` and `library_type`; the others describe the donor (`sex`, `family_id`, `father`, `mother`), the sample (`is_tumor`, `purity`) and the library (`library_kit`, `folder_name`), see {ref}`the overview <overview>`.
This one has the normal and the tumor of one patient:

```tsv
donor_name  sex     sample_name  is_tumor  library_name       library_type  library_kit
P001        female  P001-N1      false     P001-N1-DNA1-WES1  WES           Twist_Exome_2
P001        female  P001-T1      true      P001-T1-DNA1-WES1  WES           Twist_Exome_2
```

`snappy sheet check samplesheet.tsv` lists what the sheet holds, or every problem with its line and column.

The reads of a library are below `raw`, in a folder named like its `folder_name` (else its name) or named after it (`P001-N1-DNA1-WES1_S1_L001_R1_001.fastq.gz`), flat or in batch folders; the common namings of sequencing facilities and SRA are known.
`snappy scan-files` searches once and writes `reads.tsv`, one row per read group, which the tasks read from then on.
Other namings take `--pattern key=regex` for the path below the folder or the file name, whose `(?P<readgroup>...)` group pairs the `left` and `right` file of a read group, and another directory takes `--search-path`; the table records the command.
Run it again after adding libraries; `snappy run` stops for a library without rows and names the command.

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

`snappy check-files` lists the files that `config.yaml` names but that do not exist, the files of `reads.tsv` included, and the libraries without rows in it.
Then start with a dry-run.
It validates `config.yaml` and lists the jobs without running them.
Everything after `--` goes to Snakemake.

```shell
snappy run -- --dry-run --cores 1
```

Then run locally, or through the cluster executor of a Snakemake profile (see {ref}`usage`):

```shell
snappy run -- --cores 8
snappy run -- --profile slurm
```

The files of a task are in `tasks/<task_name>/work/`.
`tasks/<task_name>/output/` links its result files, see {ref}`usage`.

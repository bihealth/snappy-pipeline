(generic_pipeline_step)=

# Generic Pipeline Step Description

A pipeline step is a Snakemake workflow plus the Python code that wires its rules.
A **task** is a configured use of a step in the `tasks` list of the `config.yaml` of a project.
Each task reads files that other tasks or the project provide, writes its results into its `work/` directory and publishes the result files in its `output/` directory.

## Configuration

A task has a `step`, a unique `name` and a `config`.
The `config` is checked against the model of the step when `snappy run` starts.
Unknown keys are an error.

```yaml
- step: variant_calling
  name: calling
  config:
    depends_on:
      alignments: mapping
      reference: genome
    tool: mutect2
    mutect2:
      contamination: {enabled: false}
```

`depends_on`
: Names the upstream task for each input.
    The keys are named after the data (`reads`, `alignments`, `variants`, `reference`, `index`, ...).
    Each step page lists its keys and the data type and tags the upstream task must produce.

`tool`
: Selects the program the task runs.
    The configuration of the tool goes into the key of the same name (`mutect2` above).
    It must be present, even if all its options keep their defaults (`{}`).

`library_selection`
: Restricts the task to the libraries that match a pandas `DataFrame.query()` expression.

`group_by`
: With `cohort`, the task writes one output per cohort (pedigree) instead of one per library, where the step supports it.

`delete_intermediates`
: Deletes per-sample intermediate files once the later rules of the task have read them, see {ref}`usage`.

Each step page ends with the default configuration, with every option and its documentation.
`snappy task add` writes the same into `config.yaml`.

```{note}
**Debugging configuration errors.**
Relative paths are relative to the project directory, and `pathvars` can name roots for paths, see {ref}`config_paths`.
`snappy check-files` lists the files that `config.yaml` names but that do not exist.
Use absolute paths for data that should not move with the project.
`snappy run --verbose` logs the Snakemake command line.
`snappy run -- --dry-run` validates the configuration and lists the jobs without running them.
```

## File system layout

```text
project/
+-- config.yaml
+-- samplesheet.tsv
+-- tasks/
    +-- <task_name>/
    |   +-- work/
    |   +-- output/
    +-- <another_task>/
        +-- work/
        +-- output/
```

Each task has its own directory `tasks/<task_name>/`.
Tasks of the same step with different tools or settings do not collide.
Paths below `work/` and `output/` contain only what varies, such as the library name.
The tool is part of the task configuration and never part of a path.

There are no link-in jobs and no `input/` directory.
Steps read FASTQ files where they are, and existing data comes from `external_data` tasks.

### Directory `work/`

`work/` holds all results of the task, intermediate and final, and its logs.
Other tasks read their inputs from here.
Paths and file names in `work/` are **unstable**: the authors of a step can change them between minor versions.

Scatter fragments are temporary files, and `delete_intermediates` controls the per-sample intermediates, see {ref}`usage`.
A rule's log is its log file.
Only published files have a `.md5` file.

### Directory `output/`

`output/` is the public directory with **stable** paths and names.
It holds relative symlinks into `work/` that mirror the paths below `work/`.
No task reads from it.

One local job per task, `<task>_publish`, creates the links to the result files of the task.
Each file has its `.md5` file, and VCF and BED files also have their index.
The job writes `output/.manifest.tsv`, which lists the path, size and md5 file of each published file, and removes links to files that are no longer published.

(usage)=

# Usage

The `snappy` command is the interface to the pipeline.
`snappy --help` lists the subcommands.

## Subcommands

| Subcommand           | Purpose                                                                               |
| -------------------- | ------------------------------------------------------------------------------------- |
| `snappy init`        | Create a project directory with `config.yaml`, `pipeline_job.sh` and `README.md`.     |
| `snappy task add`    | Add tasks to the `config.yaml` of a project.                                          |
| `snappy task list`   | List the tasks of a project, or the available steps with `--available`.               |
| `snappy task steps`  | List the available steps.                                                             |
| `snappy check-files` | List the files that `config.yaml` names but that do not exist.                        |
| `snappy run`         | Run the pipeline with Snakemake.                                                      |
| `snappy report`      | Build the Snakemake report of the QC outputs (`report.zip` in the project directory). |
| `snappy watch`       | Open the `snkmt` TUI to monitor a running workflow.                                   |
| `snappy logs`        | Archive the logs of a run whose Snakemake process was killed.                         |
| `snappy refresh`     | Recreate `pipeline_job.sh` and ensure that `slurm_log` exists.                        |
| `snappy status`      | Check the status of a SLURM job with `sacct`.                                         |
| `snappy pull-sheet`  | Pull a sample sheet from the SODAR API.                                               |
| `snappy completion`  | Print the shell completion script for bash, zsh or fish.                              |

`snappy init` and `snappy task add` take `[task_name=]step_type` arguments, for example `mapping=ngs_mapping`.
`--config-mode` sets the style of the generated task configuration: `commented` (default), `full` or `minimal`.

## snappy run

`snappy run` reads `config.yaml` from the current directory, or from `--directory`.
It validates the whole project first, so a configuration error stops the run before any job starts.
By default, `snappy run` builds and publishes every task.

`--task <name>`
: Build and publish only the named task.
    Upstream tasks whose files are missing or outdated are built as well, but not published.

`--frozen`
: Needs `--task`.
    Only the rules of that task are loaded.
    The files of the upstream tasks are used as they are and never rebuilt.
    If one is missing, Snakemake stops with a missing-input error that names the file.

`--slurm`
: Run the jobs through SLURM.
    This adds the SLURM executor profile to the default profile.

`--verbose`
: Log more, including the Snakemake command line.

Everything after `--` goes to Snakemake:

```shell
snappy run -- --dry-run --cores 1               # list the jobs, run nothing
snappy run --task mapping -- --cores 8          # one task and its missing upstream files
snappy run --task calling --frozen -- --cores 8 # reuse the existing mapping, never rerun it
snappy run --slurm                              # all tasks through SLURM
```

## Directories of a task

Each task has its own directory `tasks/<task_name>/`.

`tasks/<task_name>/work/`
: Holds all files of the task, including logs.
    Other tasks read their inputs from here.
    Paths in `work/` can change between versions.

`tasks/<task_name>/output/`
: Is for publication only, and no task reads it.
    It holds relative symlinks into `work/` with the same paths, and `output/.manifest.tsv`.
    One local job per task, `<task>_publish`, creates them.
    The manifest lists the path, size and md5 file of each published file.
    Links to files that are no longer published are removed.

Only published files have a `.md5` file next to them.
Logs are not published.

## Report

`snappy report` builds the Snakemake report of the QC outputs, with one category per task (FastQC, Picard metrics, `bcftools stats`, cnvkit plots).
It uses the same Snakefile, directory and profile as `snappy run`, and arguments after `--` go to Snakemake.
`--output` sets the file below the project directory.
A `.zip` report (the default `report.zip`) can hold HTML directories such as the FastQC reports.

## Log archive

When a run ends, successfully or not, `snappy run` writes `logs/<run start>.tar.gz`.
It holds the Snakemake log of the run, the logs of the jobs that ran with their SLURM logs, and a copy of the `snkmt` database.
The Snakemake log names these files, so the task directories are not searched.
Dry-runs write no archive.

A killed Snakemake process (walltime, `scancel`) does not run its end hooks.
`snappy logs` then builds the archive from the newest Snakemake log, or from the log given as argument.

## Temporary files

Scatter fragments are the per-region or per-chunk pieces that a later rule of the task joins.
They are `temp()` files, and Snakemake deletes them when the joining job is done.

Per-sample intermediates are kept, unless the task sets `delete_intermediates: true`.
These are files that only later rules of the same task read, such as raw calls before filtering, coverage files, the normals of a panel and gVCFs.

These files are never deleted:

- Files in `output/`.
- Files that a task passes to other tasks.
- Project-wide resources, such as installed software and prepared references.

So `--task` and `--frozen` runs find what they need.

Once a temporary file is gone, a change to a later rule (a parameter or its code) also reruns the rules that produced the file.
For example, this reruns all scatter jobs of a caller.
With `delete_intermediates`, it reruns the calling when only the filter changed.

Snakemake options for temporary files, after `--`:

`--notemp`
: Keep all temporary files of this run, for debugging.

`--all-temp`
: Treat every output as temporary, to save space in test runs.

`--delete-temp-output`
: Delete temporary files left over by earlier runs.
    With `--dry-run`, only list them.

`--keep-incomplete`
: Keep the outputs of failed jobs, for debugging.

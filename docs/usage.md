(usage)=

# Usage

The `snappy` command is the interface to the pipeline.
`snappy --help` lists the subcommands.

## Subcommands

| Subcommand           | Purpose                                                                               |
| -------------------- | ------------------------------------------------------------------------------------- |
| `snappy init`        | Create a project directory with `config.yaml` and `README.md`.                        |
| `snappy task add`    | Add tasks to the `config.yaml` of a project.                                          |
| `snappy task list`   | List the tasks of a project, or the available steps with `--available`.               |
| `snappy task steps`  | List the available steps.                                                             |
| `snappy check-files` | List the files that `config.yaml` names but that do not exist.                        |
| `snappy run`         | Run the pipeline with Snakemake.                                                      |
| `snappy report`      | Build the Snakemake report of the QC outputs (`report.zip` in the project directory). |
| `snappy watch`       | Open the `snkmt` TUI to monitor a running workflow.                                   |
| `snappy logs`        | Archive the logs of a run whose Snakemake process was killed.                         |
| `snappy pull-sheet`  | Pull a sample sheet from the SODAR API.                                               |
| `snappy completion`  | Print the shell completion script for bash, zsh or fish.                              |

`snappy init` and `snappy task add` take `[task_name=]step_type` arguments, for example `mapping=ngs_mapping`.
`--config-mode` sets the style of the generated task configuration: `commented` (default), `full` or `minimal`.

## snappy run

`snappy run` reads `config.yaml` from the current directory, or from `--directory`.
It validates the whole project first, so a configuration error stops the run before any job starts.
By default, `snappy run` builds every task and links its result files into `output/`.

`--task <name>`
: Build only the named task, and link its result files.
    Upstream tasks whose files are missing or outdated are built as well, without linking their result files.

`--frozen`
: Needs `--task`.
    Only the rules of that task are loaded.
    The files of the upstream tasks are used as they are and never rebuilt.
    If one is missing, Snakemake stops with a missing-input error that names the file.

`--verbose`
: Log debug messages, too.

Everything after `--` goes to Snakemake:

```shell
snappy run -- --dry-run --cores 1               # list the jobs, run nothing
snappy run --task mapping -- --cores 8          # one task and its missing upstream files
snappy run --task calling --frozen -- --cores 8 # reuse the existing mapping, never rerun it
snappy run -- --profile slurm                   # all tasks through the executor of a profile
```

## Cluster execution

`snappy run` always uses the workflow profile of snappy, which deploys the software with conda and, for the tools that run from a container image, with apptainer, and sets the `snkmt` logger and the database persistence backend.
The executor and its settings belong to the system or the user: a Snakemake profile below `/etc/xdg/snakemake/` (system) or `~/.config/snakemake/` (user), selected with `--profile <name>` after `--` or with the environment variable `SNAKEMAKE_PROFILE`.
For SLURM, for example `~/.config/snakemake/slurm/config.yaml`:

```yaml
executor: slurm
jobs: 250
restart-times: 3
max-jobs-per-second: 1
max-status-checks-per-second: 10
local-cores: 1
latency-wait: 60
default-resources:
  slurm_partition: medium
  runtime: 1440
  mem_mb: 5000
slurm-logdir: slurm_log
slurm-keep-successful-logs: true
```

The Snakemake process itself runs as long as the workflow.
Start it in a terminal multiplexer such as `tmux` on the login node, or in an allocation (`salloc`), as the documentation of the SLURM executor recommends.

On a cluster, every job is a cluster job, unless its rule is local.
The task field `job_group` puts the jobs of a task into a Snakemake job group: dependent jobs of the tasks with the same `job_group` run as one cluster job, for example an annotation and the filters after it, per library.
Their resource requests add up (runtimes of jobs that run one after another, cores and memory of jobs that run side by side), so the members should share a partition.
`--group-components <group>=<n>` puts n independent groups into one cluster job, such as many small jobs per library:

```shell
snappy run -- --group-components tail=10
```

## Containers

PureCN and pVACtools run from container images, which their tool sections pin by digest (`purecn.container`, `pvactools.container`); a local `.sif` file works too.
Snakemake pulls each image once into `.snakemake/software/` of the project and runs the jobs of these tools inside it with apptainer.
A job in a container sees `$HOME`, `/tmp`, the project directory and the snappy installation.
Other directories that it reads, such as references outside the project, go into `APPTAINER_BIND`, a comma-separated list:

```shell
export APPTAINER_BIND=<reference directory>,<another directory>
```

apptainer keeps image layers in `~/.apptainer/cache`.
On clusters with a small home quota, point `APPTAINER_CACHEDIR` and `APPTAINER_TMPDIR` to a scratch directory.

## Directories of a task

Each task has its own directory `tasks/<task_name>/`.

`tasks/<task_name>/work/`
: Holds all files of the task, including logs.
    Other tasks read their inputs from here.
    Paths in `work/` can change between versions.

`tasks/<task_name>/output/`
: Holds only the result files, and no task reads it.
    It holds relative symlinks into `work/` with the same paths, and `output/.manifest.tsv`.
    One local job per task with result files, `<task>_results`, creates them.
    The manifest lists the path, size and md5 file of each result file.
    Links to files that are no result files any more are removed.
    A result file that is missing in `work/`, for example after a failed rerun, is built and linked again by the next `snappy run`.

Only result files have a `.md5` file next to them.
Logs are no result files.

## Report

`snappy report` builds the Snakemake report of the plots and QC outputs, with one category per task (for example FastQC, Picard metrics, BAM statistics, `bcftools stats`, the copy number plots, peddy).
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

Per-sample intermediates are kept, unless the task sets `delete_intermediates: [intermediates]`.
These are files that only later rules of the same task read, such as raw calls before filtering, coverage files, the normals of a panel and gVCFs.

These files are never deleted:

- Files in `output/`.
- Files that a task passes to other tasks.
- Project-wide resources, such as installed software and prepared references.

So `--task` and `--frozen` runs find what they need.

The one exception is the BAM of an `ngs_mapping` task that a later task processes further, such as duplicate marking or BQSR.
With `delete_intermediates: [bam]`, the BAM and its index are temporary, have no `.md5` and are not linked into `output/`; the BAM QC reports stay results.
Snakemake deletes the BAM once the tasks that read it are done, and maps again when a later run needs it.
A `--frozen` run of a later task cannot map again, so it needs the BAM to still exist.

Once a temporary file is gone, a change to a later rule (a parameter or its code) also reruns the rules that produced the file.
For example, this reruns all scatter jobs of a caller.
With `delete_intermediates: [intermediates]`, it reruns the calling when only the filter changed.

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

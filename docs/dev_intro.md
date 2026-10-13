(dev_intro)=

# Developer's Introduction

```{note}
Before reading this chapter, you should know the pipeline from the user's perspective (start at {ref}`usage`).

After reading this chapter, you should

- know which Python techniques a snappy developer needs,
- know the components of a pipeline step and how the orchestrator wires them,
- know that `snappy` is the main CLI entry point.
```

This part of the documentation targets developers who change or extend the pipeline.
It describes the architecture and points to the code of existing steps.
Most of the system follows one programming and architecture style.
New code should follow it, and {ref}`dev_conventions` states its rules.

## Prerequisites

snappy is implemented in Python (3.12 or newer) and built on [Snakemake](https://snakemake.readthedocs.io) (9 or newer).
For distributed execution, it is tailored towards SLURM.
You should be familiar with all three.

These techniques come up in the code:

- Snakemake as a rule-based language for describing workflows, including `unpack()`, scatter-gather and module namespacing,
- object-oriented programming, and classes as objects,
- Pydantic models for configuration,
- decorators, generators and closures,
- `itertools` and the built-ins `zip` and `map`,
- `snakemake.io.expand()`,
- text manipulation with `str.format`.

Useful references:

- [The Hitchhiker's Guide to Python](https://docs.python-guide.org/)
- Slatkin, Brett: Effective Python: 125 Specific Ways to Write Better Python

## Architecture

The project directory holds a `config.yaml` with a list of `tasks` and the `data_sets`.
Every task names a `step` and has a `config` for it.
`depends_on` inside a task's `config` links tasks: it names the task that provides each input.

`snappy run` calls Snakemake with the `snappy_pipeline/Snakefile`.
That file does three things:

1. It calls `load_project()` (`snappy_pipeline/orchestration.py`) once.
   This validates the config, resolves `depends_on`, rejects cycles and checks the data contracts between tasks (see {ref}`dev_contracts`).
2. It creates one workflow object per task.
3. It loads the Snakefile of each task's step as a Snakemake module.
   The module writes below `tasks/<task_name>/` (`output/`, `work/`, `benchmarks/`).

The step Snakefiles do not set up anything themselves.
They fetch their workflow object with `task_instance(config["__task_name__"])`.
By default, `snappy run` builds and publishes every task.
`--task <name>` builds and publishes one task, and `--frozen` with `--task` reads the outputs of the upstream tasks as existing files.

`work/` is the contract between tasks.
A task reads the `work/` files of its upstream tasks and never their `output/` files.
`output/` is for publication only: one local `<task>_publish` job per task links the result files into `output/` and writes `output/.manifest.tsv`.

## Anatomy of a step

For each step, there is a directory below `snappy_pipeline/workflows/` with these files:

- `__init__.py` holds the module docstring (the step's documentation page), the workflow class and the step part classes.
- `model.py` holds the Pydantic config model, including the `depends_on` model.
- `Snakefile` holds the rule definitions, which only call the workflow object.

The workflow class derives from {class}`BaseStep <snappy_pipeline.workflows.abstract.BaseStep>`.
It sets `name`, `config_model_class` and `produces`, and it registers the step parts with `register_sub_step_classes(...)`.
The workflow classes are listed in `WORKFLOW_REGISTRY` in `snappy_pipeline/workflow_registry.py`.
A new step needs an entry there.

A step part derives from {class}`BaseStepPart <snappy_pipeline.workflows.abstract.BaseStepPart>`.
It implements one tool or one helper of the step.
A step part has a tuple of `actions`.
An action is one rule of the step part.
A tool that runs as a single rule has the one action `run`.
A tool that needs several rules has several actions.
For example, mutect2 has the actions `scatter`, `run`, `gather`, `filter`, `contamination`, `pileup_normal` and `pileup_tumor`.

A step has a singular `tool` in its task config.
The workflow class registers only the step part of the selected tool, plus the structural helper step parts.
The Snakefile defines only the rules of the selected tool.
Combining several tools means defining several tasks, not multiplexing within one task.

For each action, the step part implements methods that the Snakefile reaches through the workflow object:

```python
rule ngs_mapping_bwa_run:
    input:
        unpack(wf.get_input_files("bwa", "run")),
    output:
        **wf.get_output_files("bwa", "run"),
    log:
        **wf.get_log_file("bwa", "run"),
    benchmark:
        wf.get_benchmark_file("bwa", "run")
    threads: wf.get_resource("bwa", "run", "threads")
    resources:
        runtime=wf.get_resource("bwa", "run", "runtime"),
        mem=wf.get_resource("bwa", "run", "mem"),
        partition=wf.get_resource("bwa", "run", "partition"),
        tmpdir=wf.get_resource("bwa", "run", "tmpdir"),
    params:
        args=wf.get_params("bwa", "run"),
    wrapper:
        wf.wrapper_path("bwa")
```

{ref}`dev_conventions` lists the methods and the per-action methods behind them.

Business logic belongs in the Python classes.
Python has tooling for tests, static analysis, documentation and style, and these tools stop at the first `rule` keyword of a Snakefile.
The Snakefile only wires the Python code into rules.

## Wrappers

The tools run in wrappers below `snappy_wrappers/wrappers/<name>/`.
A wrapper directory holds `wrapper.py`, `environment.yaml` and its pin file `environment.linux-64.pin.txt`.
`BaseStep.wrapper_path(...)` turns the name into the file URL for the `wrapper:` section.

A wrapper derives from a base class in `snappy_wrappers/snappy_wrapper.py`:

`ShellWrapper`
: Runs a Bash script.

`RWrapper`
: Runs an R script.

`PythonWrapper`
: Runs a Python function.

All three write the log, and they write the `<file>.md5` of every output that declares one.
A wrapper reads the values from `snakemake.params.args` and the files from `snakemake.input`.
It never reads `snakemake.config`.

## Anatomy of the `snappy` command

`snappy` is a Click application (`snappy_pipeline/apps/snappy_cli.py`).
It wraps Snakemake and reads the project configuration from `config.yaml` in the project directory.

`snappy init`
: Creates a project directory.

`snappy task add`
: Adds a task of a step to the config with the default values of its step.

`snappy task list`, `snappy task steps`
: List the tasks of the project and the available steps.

`snappy run`
: Calls Snakemake with defaults for local execution, or with `--slurm` for SLURM.
    Other arguments go to Snakemake.

`snappy watch`, `snappy status`, `snappy logs`, `snappy report`
: Monitor a run, check a SLURM job, archive the logs of a run and build the Snakemake report.

`snappy run` prints the Snakemake command line at the start of every run.
The generally useful `snakemake` parameters are available through `snappy run`.
If one is missing, please create a ticket.

(dev_intro)=

# Developer's Introduction

```{note}
Before reading this chapter, you should know the pipeline from the user's perspective (start at {ref}`usage`).

After reading this chapter, you should

- know which Python techniques a snappy developer needs,

- know the components of a pipeline step and how the orchestrator wires them,

- know that `snappy` is the main CLI entry point.
```

This part of the documentation is for developers who change or extend the pipeline.
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

### From `config.yaml` to rules

The project directory holds a `config.yaml` with a list of `tasks` and the `sample_sheets`.
Every task names a `step` and has a `config` for it.
`depends_on` inside a task's `config` links tasks: it names the task that provides each input.

`snappy run` calls Snakemake with the `snappy_pipeline/Snakefile`.
That file does three things:

1. It calls `load_project()` (`snappy_pipeline/orchestration.py`) once.
   This validates the config, resolves `depends_on`, rejects cycles and checks the data contracts between tasks (see {ref}`dev_contracts`).
2. It creates one workflow object per task.
3. It loads the Snakefile of each task's step as a Snakemake module.
   The module writes below `tasks/<task_name>/` (`output/`, `work/`, `benchmarks/`).

```{graphviz}
:alt: config.yaml, load_project(), one workflow object per task, one Snakemake module per task, and the generated and explicit rules

digraph config_to_rules {
    rankdir=LR; newrank=true; nodesep=0.12; ranksep=0.45;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    config [label=<<B>config.yaml</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tasks and sample_sheets</FONT><BR ALIGN="LEFT"/>>];
    subgraph cluster_snakefile {
        label=<<FONT POINT-SIZE="7" COLOR="#898781">snappy_pipeline/Snakefile</FONT>>; color="#e1e0d9"; style="rounded"; labelloc=t;
        load [label=<<B>load_project()</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">validates the config, resolves depends_on,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">checks the data contracts</FONT><BR ALIGN="LEFT"/>>];
        objects [label=<<B>workflow objects</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">one per task</FONT><BR ALIGN="LEFT"/>>];
        modules [label=<<B>Snakemake modules</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the step Snakefile of each task,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">writing below tasks/&lt;task_name&gt;/</FONT><BR ALIGN="LEFT"/>>];
    }
    generated [label=<<B>generated rules</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">wrapper_rules.rules: one per action</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">with a wrapper</FONT><BR ALIGN="LEFT"/>>];
    explicit [label=<<B>explicit rules</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">containers, checkpoints,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">run: and shell: blocks</FONT><BR ALIGN="LEFT"/>>];
    config -> load -> objects;
    objects -> modules [label="task_instance()"];
    modules -> generated;
    modules -> explicit;
}
```

A step Snakefile neither loads the config nor builds its workflow object.
It fetches the object with `task_instance(config["__task_name__"])`.

### `work/` and `output/`

By default, `snappy run` builds every task and links its result files into `output/`.
`--task <name>` does so for one task.
With `--task`, `--frozen` reads the outputs of the upstream tasks as existing files and never rebuilds them.

`work/` is the contract between tasks.
A task reads the `work/` files of its upstream tasks and never their `output/` files.
`output/` holds only the result files.
Each task with result files has one local `<task>_results` job, which links them into `output/` and writes `output/.manifest.tsv`.

## Anatomy of a step

Each step has a directory below `snappy_pipeline/workflows/` with these files:

- `__init__.py` holds the module docstring (the step's documentation page), the workflow class and the step part classes.
- `model.py` holds the Pydantic config model, including the `depends_on` model.
- `Snakefile` holds the rule definitions, which only call the workflow object.

```{graphviz}
:alt: The config model, the workflow class, its step parts, their actions and one rule per action

digraph step_anatomy {
    rankdir=TB; nodesep=0.3; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    model [label=<<B>config model</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">model.py, with the depends_on model</FONT><BR ALIGN="LEFT"/>>];
    workflow [label=<<B>workflow class</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">__init__.py, derives from BaseStep</FONT><BR ALIGN="LEFT"/>>];
    parts [label=<<B>step parts</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">__init__.py, derive from BaseStepPart</FONT><BR ALIGN="LEFT"/>>];
    actions [label=<<B>actions</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">such as scatter, run, gather</FONT><BR ALIGN="LEFT"/>>];
    rules [label=<<B>rules</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">wrapper_rules.rules and the Snakefile</FONT><BR ALIGN="LEFT"/>>];
    { rank=same; model; workflow; }
    model -> workflow [label="config_model_class"];
    workflow -> parts [label="tool_parts[tool], helper_parts"];
    parts -> actions;
    actions -> rules [label="one rule each"];
}
```

### Workflow class and step parts

The workflow class derives from {class}`BaseStep <snappy_pipeline.workflows.abstract.BaseStep>`.
It sets `name`, `config_model_class` and `produces`.
It names its step parts in `tool_parts` and `helper_parts`, and `BaseStep.__init__` registers them.
A step whose parts depend on more than the tool registers them itself with `register_sub_step_classes(...)`.
`WORKFLOW_CLASSES` in `snappy_pipeline/workflow_registry.py` maps each step name to the name of its workflow class.
A new step needs an entry there; `WORKFLOW_REGISTRY` imports the class from `snappy_pipeline/workflows/<step name>/` on first access.

A step part derives from {class}`BaseStepPart <snappy_pipeline.workflows.abstract.BaseStepPart>`.
It implements one tool or one helper of the step.

A step part has a tuple of `actions`.
An action is one rule of the step part.
A tool that runs as a single rule has the one action `run`.
A tool that needs several rules has several actions.
For example, mutect2 has the actions `scatter`, `run`, `gather`, `filter`, `contamination`, `pileup_normal` and `pileup_tumor`.

The task config of a step names a single `tool`.
The workflow class registers only the step part of the selected tool (`tool_parts`) and the structural helper step parts (`helper_parts`).
So only that tool's rules exist.
To combine several tools, define one task per tool.

### Rules

For each action, the step part implements methods that its rule reaches through the workflow object.
The rule of every action with a wrapper comes from the shared include `snappy_pipeline/workflows/common/wrapper_rules.rules`, not from the step Snakefile:

```python
def wrapper_rule(spec):

    rule:
        name:
            spec.name
        input:
            unpack(wf.get_input_files(spec.part, spec.action)),
        output:
            **wf.get_output_files(spec.part, spec.action),
        log:
            **wf.get_log_file(spec.part, spec.action),
        benchmark:
            wf.get_benchmark_file(spec.part, spec.action)
        threads: wf.get_resource(spec.part, spec.action, "threads")
        resources:
            **wf.get_resources(spec.part, spec.action),
        params:
            args=wf.get_params(spec.part, spec.action),
        wrapper:
            wf.get_wrapper(spec.part, spec.action)


for spec in wf.wrapper_rules():
    wrapper_rule(spec)
```

The rule sits in a function because Snakemake evaluates the `wrapper:` expression again when a job runs.
A rule defined in the loop itself would run the wrapper of the loop's last value.

A step Snakefile defines `wf`, includes that file and writes only the rules the include cannot express: containers, checkpoints, `run:` and `shell:` blocks.
{ref}`dev_conventions` lists the methods and the per-action methods behind them.
It also lists the `rule_options` of a part (`localrule`, `retries`, `cache`, `wildcard_constraints`, `explicit`).

Business logic belongs in the Python classes.
Python has tooling for tests, static analysis, documentation and style.
These tools stop at the first `rule` keyword of a Snakefile.

## Wrappers

The tools run in wrappers below `snappy_wrappers/wrappers/<name>/`.
A wrapper directory holds `wrapper.py`, `environment.yaml` and its pin file `environment.linux-64.pin.txt`.
A step part names the wrapper of each action in `wrappers`.
`wf.get_wrapper(...)` turns the name into the file URL for the `wrapper:` section.

```{graphviz}
:alt: The step part names the wrapper, the rule passes inputs and params to wrapper.py, which writes the outputs, the log and the md5 files

digraph wrapper_job {
    rankdir=LR; nodesep=0.12; ranksep=0.45;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    part [label=<<B>step part</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">wrappers: the wrapper of each action</FONT><BR ALIGN="LEFT"/>>];
    rule [label=<<B>rule</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">input:, params: args=..., wrapper:</FONT><BR ALIGN="LEFT"/>>];
    wrapper [label=<<B>wrapper.py</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">ShellWrapper or PythonWrapper</FONT><BR ALIGN="LEFT"/>>];
    outputs [label=<<B>outputs</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the tool's files, the log</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">and &lt;file&gt;.md5</FONT><BR ALIGN="LEFT"/>>];
    part -> rule [label="wf.get_wrapper():\nfile:// URL"];
    rule -> wrapper [label="snakemake.input,\nsnakemake.params.args"];
    wrapper -> outputs;
}
```

A wrapper derives from a base class in `snappy_wrappers/snappy_wrapper.py`:

`ShellWrapper`
: Runs a Bash script.
    R code runs from a heredoc in the script.

`PythonWrapper`
: Runs a Python function.

Both write the log and the `<file>.md5` of every output that declares one.
A wrapper reads the values from `snakemake.params.args` and the files from `snakemake.input`.
It never reads `snakemake.config`.

## Anatomy of the `snappy` command

`snappy` is a Click application (`snappy_pipeline/apps/snappy_cli.py`).
It wraps Snakemake and reads the project configuration from `config.yaml` in the project directory.
{ref}`usage_commands` lists its commands.
`docs/conf.py` generates that table from the help of each command.

`snappy run` calls Snakemake with snappy's workflow profile.
The executor, such as SLURM, comes from a profile of the system or the user.
Arguments that `snappy run` does not define go to Snakemake.

`snappy run --verbose` prints the Snakemake command line at the start of the run.
Every argument after `--` goes to Snakemake, so its options are available through `snappy run`.
If one is missing, create a ticket.

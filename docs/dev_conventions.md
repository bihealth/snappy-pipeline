(dev_conventions)=

# Developer Conventions

These rules apply to every workflow step.
Tests enforce most of them.
The section on each rule names the check.
`tests/snappy_pipeline/test_conventions.py` holds the static checks.

## Configuration is explicit

A configuration must mean the same thing across snappy versions and setups.
Code does not guess.
There is no automatic wiring of `depends_on` and no fallback that picks a task by its step type.
Validation fails with a clear message instead of choosing for the user.

- `tool` has no default.
  Every task names its tool (`test_tool_is_required`).
- `depends_on` keys have no task-name defaults (see {ref}`Contracts <dev_contracts>`).
- A field's description (its attribute docstring) is plain text with literals in double backticks.
  It has no Sphinx roles, directives or references to names in the code.
  It becomes a YAML comment in the configuration that `snappy task add` writes.
  There, `default_config_yaml_string` wraps it and turns the double backticks into single ones.
  It also becomes a cell of the inputs table of the step's docs page (`test_field_descriptions_are_plain_text`).
- In the docs and in the test pipelines, a whole `config.yaml` lists `pathvars`, `data_sets` and `resources` before `tasks`, like `snappy init` writes it.
  A task added later then ends the file (`test_data_sets_precede_the_tasks`).

## Step part API

The orchestrator loads the project once (`snappy_pipeline.orchestration.load_project`) and creates the workflow object of every task once.
A step Snakefile fetches its workflow object instead of building one:

```python
from snappy_pipeline.orchestration import task_instance

wf = task_instance(config["__task_name__"])
```

`wf.get_task_config(name)` returns the task's own config (`name` is its step or task name).
It also returns the config of the upstream task set in `depends_on.<name>`.
Anything else is an error.

### Generated rules

A step Snakefile defines `wf` and its global `wildcard_constraints`.
Then it includes `snappy_pipeline/workflows/common/wrapper_rules.rules`.
For every registered step part, the include writes one rule per action that has a wrapper (`BaseStep.wrapper_rules`).
The rule is named `<step>_<part>_<action>`.
Snakefiles stay thin and call the workflow object for everything.

A part's `rule_options` give an action's `localrule`, `retries`, `cache` and `wildcard_constraints`.
A rule with `cache` has no benchmark.
`rule_options` can also mark an action `explicit`.
Then the Snakefile writes the rule itself.
It does so for containers (`software:`), checkpoints, `run:` and `shell:` blocks and outputs that are no dict.

### Methods

A step part provides these methods.
The generated and the explicit rules reach them through the workflow object (`wf`).
Each `wf` method takes the step part name before the action, for example `wf.get_input_files("bwa", "run")`.

| Method                       | Returns                                                                    | Rule section                                                                            |
| ---------------------------- | -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `get_input_files(action)`    | a function of wildcards, returning a dict (named inputs), a list or a path | `input: unpack(wf.get_input_files(...))` for dicts, otherwise `wf.get_input_files(...)` |
| `get_output_files(action)`   | static paths (Snakemake does not accept functions here)                    | `output: **wf.get_output_files(...)`                                                    |
| `get_log_file(action)`       | static paths                                                               | `log: **wf.get_log_file(...)`                                                           |
| `get_params(action)`         | a function of wildcards (optionally also of `input`), returning a dict     | `params: args=wf.get_params(...)`                                                       |
| `get_resource(action, name)` | a function of wildcards, input, threads and attempt                        | `threads: wf.get_resource(..., "threads")`                                              |
| `get_resources(action)`      | the functions of runtime, mem, partition and tmpdir                        | `resources: **wf.get_resources(...)`                                                    |
| `get_wrapper(action)`        | the `file://` URL of the action's wrapper, from the part's `wrappers`      | `wrapper: wf.get_wrapper(...)`                                                          |

`BaseStep` raises a `TypeError` at parse time if `get_input_files` or `get_params` returns anything but a function.

A step part implements these methods through per-action methods.
`BaseStepPart` validates the action and dispatches to them:

- `_get_input_files_<action>(wildcards)`
- `_get_params_<action>(wildcards)`
- `_get_output_files_<action>()`
- `_get_log_file_<action>()`

One `_get_log_file(action)` can replace the log methods of all actions.

### Defining input files and params

`BaseStepPart.get_input_files(action)` and `get_params(action)` validate the action against the class's `actions` tuple.
They return `self._get_input_files_<action>` or `self._get_params_<action>`.
A step part therefore only defines these methods:

```python
class ExampleStepPart(BaseStepPart):
    name = "example"
    actions = ("run",)

    def _get_input_files_run(self, wildcards):
        return {"vcf": "work/{library_name}/out/{library_name}.vcf.gz".format(**wildcards)}

    def _get_params_run(self, wildcards):
        return {"min_af": self.config.example.min_af}
```

Rules:

- `actions` is a tuple.
- Input and params methods take exactly `(self, wildcards)`.
  Params may also take `(self, wildcards, input)`.
  No `**kwargs`: Snakemake passes every available argument to a function that accepts `**kwargs`.
- Override `get_input_files` or `get_params` only if the override does more than dispatching.
  It must still return a function.
- `@dictify` and `@listify` keep the decorated function's signature (`functools.wraps`), which Snakemake inspects.
- Paths returned by input functions are relative (`work/...`, or upstream paths from `get_upstream_paths()`).
  Snakemake adds the task prefix to function results as it does to static patterns.

The exception to these signatures is an input function that needs Snakemake's `checkpoints` object.
The Snakefile calls it as `wf.get_input_files(...)(wildcards, checkpoints)` from a Snakefile-level function (`gcnv_model`).
Each task is a Snakemake module with its own `checkpoints` object, so a task cannot use the checkpoint of an upstream task.
It declares a local checkpoint instead (`sv_calling_*` with `tool: gcnv`).
The local checkpoint reads the output of the upstream checkpoint and lists what the upstream checkpoint produced.

### Registered parts and wrappers

A step names the part of each tool in `tool_parts` and its structural parts (`write_pedigree`) in `helper_parts`.
`BaseStep.__init__` registers the helpers and the selected tool's part.
A step registers its parts itself only if they depend on more than the tool (`register_sub_step_classes(...)`).

A step part's `wrappers` maps each action that runs a wrapper to the wrapper's path below `snappy_wrappers/wrappers`.
Where the config picks the wrapper, `wrappers` is a property.
Examples are the whole-genome variants of the gCNV wrappers and the tool of a filtration task.
Every rule that runs a wrapper passes `params: args=wf.get_params(...)`.
An action without params of its own gets an empty dict, to which `BaseStep.get_params` adds `extra_arguments` and `code`.

### Module and class layout

A step module lists its members in this order:

1. its imports, constants and functions,
2. the step parts in the order the step class registers them (a base part before its subclasses),
3. the step class.

A step or step part class lists its attributes first.
The attributes of its base classes come first, in the bases' order, root first.
So `name`, `actions`, `wrappers`, ... come before a part's own attributes.
A property that overrides a base attribute counts as that attribute.
Then come the methods:

1. `__init__` and the other overrides of base methods, in the bases' order,
2. the per-action methods by action, in `actions` order, and within an action inputs, outputs, log file and params,
3. the class's own methods.

`pixi run lint-layout` checks this, and `pixi run layout` rewrites the modules (`tests/scripts/class_layout.py`).
`test_workflow_modules_follow_the_layout` fails otherwise.

### Resource usage and result files

Resource usage is declared, not computed.
`default_resource_usage` holds it for all actions, and `resource_usage` overrides it per action.
Either is a property when it depends on the config.

`result_keys` names the output keys whose files are result files, for all actions or by action.
`BaseStepPart.get_result_files` derives their `work/` paths for each output entity (`entity_wildcard`, other wildcards from `_result_wildcards`).
It derives a static path, such as a panel of normals, once.
`BaseStep.get_result_files` collects them from all step parts.
Override these only where the derivation does not fit.

### Libraries, donors and pedigrees

Step parts read the sample sheets only through the library DataFrame of the workflow object, never through the biomedsheets objects.
`build_library_dataframe()` returns one row per library.
A row holds the library's role, donor, group of related libraries (`related_group`), parents, sex and the extra sample sheet columns.
`include_background=True` adds the libraries of background data sets.
The project loads each sample sheet once (`Project.data_set_infos`).
The table without relationships is built once per run and shared by all tasks (`BaseStep.library_table()`, `Project.library_tables`).
Each task adds its `relationships` to it.

For germline and generic sheets, the table lists all libraries with the pedigree columns of the biomedsheets germline view.
That view joins the pedigrees by the data set's `pedigree_field`, if set.
`library_row(name)` returns the row of one library.
If several data sets have a library of that name, `library_selection` picks the row.
Read per-library facts (data set, folder, kit) through it, not by indexing the dataframe by name.

`libraries` is the one selection of the libraries a task takes:

1. the library table, with background data sets if the step sets `background_results`,
2. restricted to the step's `library_requirement`,
3. filtered by `library_selection`, or else by the step's `library_default`.
   A `library_default` of `None` takes the default of each library's sheet kind (cancer: tumor DNA, germline: DNA), and `""` takes all libraries.

`library_requirement` is what a step can process at all (RNA for expression quantification, `kind == 'cancer'` for the somatic steps).
It always applies.
`library_selection` is the user's choice within it, such as a batch or dropped samples.
`library_default` replaces it when the task sets none.
A step whose requirement or default depends on the task's config (its tool) overrides `task_library_requirement(config)` or `task_library_default(config)`.

Steps and tools that compare each tumor with its matched normal (the `matched_normal_lib` relationship) require `MATCHED_NORMAL_REQUIREMENT`.
They leave out tumors without one, and `BaseStep.library_warnings()` names these libraries.
Tumor-only tools keep them.

Only the default table of a step that sets `background_results` has the background libraries.
Steps that set `uses_background` add them where they need them.
A task's `background_selection` narrows them (`build_library_dataframe(include_background=True)` and `task_libraries` apply it).

`task_libraries(task_name)` gives the libraries of another task, computed from its config.
So a step can match the files of an upstream task even under `--frozen`, which creates no workflow objects for upstream tasks.
`output_entities` gives the entities a task writes, from `libraries` and `effective_group_by`.
`effective_group_by` is the task's `group_by`, else the step's `default_group_by(config)`, else what upstream tasks declare, else `library`.
`primary_library(entity)` gives the library that stands for an entity: itself, or the `is_primary` member of its related group.
`member_libraries(entity)` gives the libraries that an entity stands for in a joint result.
These are the pedigree members, one DNA library per donor, and only an explicit `library_selection` restricts them.
A step never treats an entity as a library.

`primary_dna_libraries()`, `pedigrees()`, `pedigree_members()`, `pedigree_of(library)`, `pedigree_index(pedigree)` and `children_with_both_parents()` give the germline view: one DNA library per donor, grouped by pedigree.
`pedigree_donors()` gives the donor table of a pedigree (`build_donor_dataframe()`), for steps that model unsequenced relatives.
The donor table has every donor of the germline sheets, sequenced or not.
A row holds the donor's parents as donor names, its sex and its DNA library.

## Wrappers

### Params

Wrappers read the params dict as `snakemake.params.args`.
Snakemake does not allow `unpack()` in `params:`, so every rule passes it as the single key `args`.

Wrappers never read `snakemake.config`.
Every file a wrapper reads, such as the reference FASTA, is a rule input, and every other value is a param.
Only inputs let Snakemake order the jobs and rerun them when a file changes.
The job that produces an input runs first, and a changed input reruns the job that reads it.
A changed param also triggers a rerun.
So every value a wrapper uses is a config field passed as a param, never a constant in the wrapper or the Snakefile.
Rules with a `run:` block need params for the config values their output depends on, too.

- `test_wrappers_get_the_params_they_read` checks that the jobs in the DAG snapshots pass every `args` key their wrapper reads.
- `test_wrappers_read_the_params_they_get` checks the other direction: every param of a wrapper job occurs in its wrapper.
  `WRAPPERS_READING_PARAMS_DYNAMICALLY` lists the exceptions with a reason.
- `test_files_are_inputs_not_params` fails for params that hold file paths.

### Code hash

The params of every rule that runs a wrapper carry `code`, a hash of the wrapper's code (`snappy_pipeline/wrapper_code.py`).
The hash covers:

- the files of the wrapper's directory,
- the files they import from `snappy_wrappers`,
- the scripts of the package they run by path or as a console script.

`wrapper_code.py` finds the files by reading the imports, without importing the tools.
Python files are hashed without docstrings and comments, so an edit to those reruns nothing.
Snakemake's own code trigger covers shell commands and `run:` blocks only and records nothing for a `wrapper:` rule.
So the `code` param reruns the jobs of a changed wrapper, tool or base class (`test_wrapper_jobs_carry_the_code_of_their_wrapper`).

### Extra arguments

The config section of a tool inherits `ToolModel` (`snappy_pipeline/models`).
Its `extra_arguments` are additional arguments for the tool's command line.
`BaseStep.get_params` takes them from the step part's `tool_section()` (by default `config.<name>`).
It adds them to the params of the part's `extra_arguments_action`, or else of its only action.
If the params method returns `extra_arguments` itself, its value wins.
The wrapper puts `{extra_arguments}` on the tool's command, and `ShellWrapper` fills the placeholder in.
A wrapper that builds its command in an f-string writes `{{extra_arguments}}`, which leaves the placeholder to the base class.

A single-tool step without a tool section takes `extra_arguments` on its step model.
Its step part's `tool_section()` returns `self.config`.
`test_tool_wrappers_get_extra_arguments` requires them for every wrapper.
`WRAPPERS_WITHOUT_EXTRA_ARGUMENTS` lists the exceptions with a reason: Python and R code, installers, helpers and secondary actions.
`test_wrappers_use_the_extra_arguments_they_get` and `test_wrappers_leave_extra_arguments_to_the_base_class` check that wrappers use them.

### Shared shell text

Several wrappers share shell text, such as the GNU parallel sharding of the GATK callers, the bwa-mem pipeline and the JVM heap.
The shared text comes from functions in `snappy_wrappers/tools`.
Some of them return finished shell text, which the wrapper pastes into its template.
Others run the whole wrapper (`bwa_mem.run`).

The convention tests read a wrapper's own files and the tool modules behind it (`wrapper_code.code_files`).
`tests/scripts/shell_templates.py` checks the shared shell text with shfmt as well.
It also flags a `set -x` or `set -euo pipefail` on an unindented line of a template, which the `ShellWrapper` header already sets.
The `ShellWrapper` keeps the script it ran next to the job log, `<name>.sh` for `<name>.log`, so `bash` reruns the job.
A wrapper embeds no script in a heredoc: `ShellWrapper.run(cmd, R=code)` writes the R code next to the log as `<name>.R` with its placeholders filled in and names it in `cmd` as `{scripts.R}` (`py=code` and `{scripts.py}` for Python), and `RWrapper.run(code)` runs an R script that way with `Rscript --vanilla`.
So the R or Python code of a job is a file to read and to run on its own.

### Environments

Every wrapper environment is pinned.
Snakemake builds it from the exact package URLs in `environment.linux-64.pin.txt` next to `environment.yaml`.
After changing an environment file, pin it again with `pixi run -e pins update-pins <environment.yaml> ...` (snakedeploy, needs conda).
Without file arguments, the task pins every environment that has a pin file.
A symlinked `environment.yaml` gets a symlink to its target's pin file.

The `conda-autobump` GitHub workflow updates and re-pins the environments weekly, with one pull request per environment.
`test_wrapper_envs_are_pinned` checks that each environment has its pin file.
It also checks that the pin file holds every package of the environment file, at the exact versions the environment file names.
In a pull request, CI also checks that a changed environment file comes with its changed pin file.

Some tools would download an asset in every job, such as igv.js for igv-reports.
The environment's `environment.post-deploy.sh` fetches such an asset once into `$CONDA_PREFIX/share/<tool>/` and checks its checksum.
The wrapper reads it from `sys.prefix`.
The calls report inlines it into the IGV page (`snappy_wrappers/tools/calls_report/igv.py`).

## Outputs

### Benchmarks

Rules that run a wrapper declare `benchmark: wf.get_benchmark_file(part, action)`.
The benchmark file has the path of the rule's log, below `benchmarks/` (`tasks/<task>/benchmarks/` in a project).
Snakemake writes the wall time, memory and I/O of each job there.
These files calibrate the resource requests.

Rules cached between workflows (`cache:`) cannot have a benchmark.
`test_wrapper_rules_write_a_benchmark` and `test_logs_and_benchmarks_are_unique_per_job` check this.

### Logs

A rule's log holds only its log file (`test_logs_hold_only_logs`).
The log file is named `log`, because the wrapper base classes write to `snakemake.log.log`.
Building the DAG in the snapshot tests fails for a wrapper rule without a log named `log`.

Snakemake's metadata in `.snakemake/` records the code, params, inputs and conda environment of every output.
The pinned environments (`environment.linux-64.pin.txt`) fix the package versions.
Logs are no result files.
`snappy run` archives them per run (`test_logs_are_no_result_files`).

### Result files and `output/`

Result files are the files linked into `output/`.
Only result files have an md5 file, and every result file has one, except directories and the results of the steps in `RESULTS_WITHOUT_MD5` (`reference_index`, `gcnv_model`, `cbioportal_export`).
SODAR needs an md5 file per uploaded file.
A step declares `<file>.md5` next to `<file>` in the outputs, and the wrapper base classes write it.
Wrappers never read md5 files.

A step's `get_result_files` lists the result files as `work/` paths, each with its md5 file.
Bgzipped VCF and BED files also come with `.tbi` and `.tbi.md5`.
Each task with result files has one local job, `<task>_results` (`snappy_pipeline/results.py`).
It links the result files and writes `output/.manifest.tsv`.
Only the results job derives their `output/` links, and no other rule writes into `output/`.

```{graphviz}
:alt: The results job links the result files from work/ into output/; downstream tasks read work/; rule all targets the result files and the manifest

digraph outputs {
    rankdir=LR; newrank=true; nodesep=0.2; ranksep=0.45;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    subgraph cluster_task {
        label=<<FONT POINT-SIZE="7" COLOR="#898781">tasks/&lt;task&gt;/</FONT>>; color="#e1e0d9"; style="rounded"; labelloc=t;
        work [label=<<B>work/</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the outputs of the jobs,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">result files with their .md5</FONT><BR ALIGN="LEFT"/>>];
        results [label=<<B>&lt;task&gt;_results</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">one local job per task</FONT><BR ALIGN="LEFT"/>>];
        output [label=<<B>output/</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">links to the result files,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">.manifest.tsv</FONT><BR ALIGN="LEFT"/>>];
    }
    downstream [label=<<B>downstream task</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">get_upstream_paths()</FONT><BR ALIGN="LEFT"/>>];
    all [label=<<B>rule all</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">targets</FONT><BR ALIGN="LEFT"/>>];
    { rank=same; results; downstream; }
    { rank=same; output; all; }
    work -> results [label="result files"];
    results -> output;
    work -> downstream [label="get_output_paths()"];
    work -> all [label="result files"];
    output -> all [label=".manifest.tsv"];
}
```

The result files in `work/` are targets of `rule all` as well.
If an input of a job is missing but the job's output exists, Snakemake does not rerun the job.
So a manifest as the only target would hide result files that a failed rerun deleted (`test_result_files_are_targets`).
`get_output_paths` returns `work/` paths, so tasks read each other's `work/` files and never `output/`.

`test_only_the_results_job_writes_into_output`, `test_no_job_reads_from_output`, `test_md5_files_only_for_result_files`, `test_result_files_come_with_md5_and_index` and `test_wrappers_read_no_md5_or_conda_files` check this.

### Directory outputs

Some tools write a directory of files that the step does not name.
Such a tool has a `directory()` output, not a `.done` marker file in the directory.
A job whose files the step names declares them and no marker next to them.
Snakemake removes the directory before the job runs and after a failure, so no file of an earlier run survives.
Snakemake creates only the parent directory, so the wrapper creates the directory itself.
No other output can be inside it, and only one job writes it.

Consumers take the directory as their input.
After the job, the directory also holds Snakemake's hidden `.snakemake_timestamp` file.
A directory can be a result file: the results job links the directory, which has no md5 file.
The DAG snapshots list the directory outputs of a job as `directory`.

### Temp files and intermediates

A step part lists output keys in these attributes.
`BaseStep.get_output_files` marks the outputs `temp()` as the table shows:

| Outputs                                                                                                                           | Step part attribute        | `temp()`                                                                                |
| --------------------------------------------------------------------------------------------------------------------------------- | -------------------------- | --------------------------------------------------------------------------------------- |
| Scatter fragments                                                                                                                 | `temp_outputs`, per action | always                                                                                  |
| Per-sample intermediates that only later rules of the same task read                                                              | `intermediate_outputs`     | when the task's `delete_intermediates` names `intermediates`                            |
| Result files, contract outputs, files other tasks read and project-wide resources (installs, prepared references, pedigree files) | neither                    | no                                                                                      |
| The exception: a BAM that another task processes further                                                                          | `bam_outputs`              | when `delete_intermediates` names `bam`, and then it is no result file (without `.md5`) |

Only steps with `bam_intermediate = True` in their config model accept `bam`.
Every temp file needs a declared reader.
Otherwise, Snakemake deletes it right after its job.
If a tool reads the index of a temp file implicitly, its job still declares the index as input.

The DAG snapshots are built with `delete_intermediates: [intermediates]`.
`test_intermediates_are_temp` and `test_temp_files_are_intermediates` check these snapshots.

### Report outputs

A step part declares the outputs for the Snakemake report per action in `report_outputs`.
It maps each output key to a `ReportOutput`: a subcategory, a caption and, for a directory, the `htmlindex` or the `patterns` of the files to show.
`BaseStep.get_output_files` flags these outputs with `report()`.
Each gets the task's category and the wildcards of its path and patterns as labels.
Captions are RST files in `report/` next to the step module.

Every plot (PDF, PNG, SVG, HTML) is a report output, and so is every `report/` directory and every output below one.
`OUTPUTS_NOT_IN_REPORT` lists the exceptions with a reason.
`test_plots_and_diagnostic_outputs_are_in_the_report` checks the DAG snapshots, which list them as `report`.

### QC summaries

The Snakemake report embeds only images, PDF and HTML.
So a step whose tool writes tables of statistics adds a QC summary.
The QC summary is a subclass of `QcSummaryStepPart` (`snappy_pipeline/workflows/common/qc_summary.py`) next to the tool part.
The subclass names:

- the `source` part and action,
- the output keys to read by kind (`source_keys`),
- the `format`, a parser module of `snappy_wrappers.tools.qc_summary` (polars).

One job per task runs the `qc_summary` wrapper.
It reads the files of all entities (`entities()`, the task's output entities by default).
`label_wildcard` labels several files of one entity, such as the samples of a VCF.
The job writes the summary table `qc_summary/out/qc_summary.tsv` (one row per library, a result file).
It also renders a datavzrd site below `qc_summary/report/html/` (a report output with `htmlindex`).

A format module turns the files into a `Summary`: the summary `Table` and the views of the site.
The views are `Table` and `Chart` objects, and a `Chart` is a vega-lite chart of a long frame, thinned to `MAX_POINTS` per library.
A check that compares the files with the sample sheet gets the sheet facts as keyword arguments of `summarize`.
They come from the part's `summary_options()` through the params, so a sheet change reruns the job.
Examples are the fingerprint checks of `ngs_chew`.
These checks get the libraries with donor, pedigree and sex, and the donors of the pedigrees with their parents.

`tests/snappy_wrappers/tools/test_qc_summary.py` checks each format against fixtures and the generated configuration against datavzrd's schema.
`test_qc_summary.py` of the workflows checks the step part.

#### Tables

A `Column` has a datavzrd display mode (`display`: `pinned`, `normal`, `available`, `detail` or `hidden`), `links` and an `ellipsis`.
`links` maps a name to a URL with `{value}` and `{Other column}` placeholders, which datavzrd fills without encoding.
A `Column` can also carry a vega-lite plot per cell (`custom_plot`: a JS `data` function of the cell value and a `spec`), or a JS renderer of the cell's HTML (`custom`), as the consequence pills of the calls report do to wrap within a bounded width.

| Function    | Shows                                                    |
| ----------- | -------------------------------------------------------- |
| `pills()`   | categories, colored on a fixed domain                    |
| `heatmap()` | a magnitude, colored on a ramp (`SEQUENTIAL`, `AGAINST`) |
| `ticks()`   | the value of a cell as a tick on the range of its column |
| `bars()`    | the value of a cell as a bar                             |

The `STATUS` colors mark states and always come with the state's label.

A `Table` can:

- fix `max_in_memory_rows`,
- link its rows to the rows of another table (`links`),
- turn off the Excel download,
- leave columns out of the site's dataset while the result TSV keeps them (`exclude`).

`Summary.default_view` names the view the site opens with.

#### Charts

A table whose values can be drawn as a chart gets one.
The table stays and shows the exact values.
The charts follow Tufte (`CHART_CONFIG`): no gridlines, hairline axes, thin lines and no frame, so the ink goes to the data.

Every summary gets a strip chart of its key metrics (`metric_strips()`).
It has one tick per library and metric, and each metric has its own scale.
So an outlier stands out, and the tooltip names it.
Per-contig counts are dot charts (`mark="dots"`).
`Chart(mark="scatter")` puts one dot per row on two quantitative axes, such as the relatedness of peddy's pairs.
Its `tooltip` names further columns for the tooltip.
`thin()` bins distributions at round widths (`nice_width()`).

Up to `MAX_COLORED` libraries get a color each, in a fixed, colorblind-safe order, with a legend that toggles them.
A distribution or dot chart of more libraries becomes a heatmap (`Chart.many`), which shows hundreds of libraries at once.
The heatmap has one row per library, the bins as columns and the value as a sequential color.

The charts are interactive:

- Continuous axes zoom and pan.
- Hovering a mark highlights its library, across the panels of a strip chart.
- A mark links to the row of its library in the summary table (`Chart.link`, a datavzrd `table-row` link named `LINK` with an `href` encoding).

### Calls report

The calls report (`snappy_wrappers/tools/calls_report/`) uses `custom_plot` for the observations and the allele frequency density of each sample.
The `calls_report` wrapper renders it per VCF for the `variant_report` step (`snappy_pipeline/workflows/variant_report/`), which reads the VCFs of `depends_on.variants`.
The columns follow the VCF's header: varlociraptor's `PROB_*` fields or the FORMAT fields of other callers.
The cell holds the compact strings of the VCF (`SAOBS|SROBS`, a thinned `AFD`).
The `data` function parses them in the browser, so a row costs tens of bytes per sample.

(dev_contracts)=

## Contracts

`DataSignature` objects describe a task's outputs: a data type plus tags such as `dna`, `somatic` or `filtered`.

```{graphviz}
:alt: load_project() checks the signatures of the upstream task against the requirement of the depends_on field; get_output_paths() and get_upstream_paths() give the consumer the work/ paths of the picked signature

digraph contract {
    rankdir=LR; nodesep=0.12; ranksep=0.45;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781", fontname="Helvetica", fontsize=7, fontcolor="#52514e"];
    upstream [label=<<B>upstream task</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">produces or</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">task_produces(config, upstream)</FONT><BR ALIGN="LEFT"/>>];
    field [label=<<B>depends_on field</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">DataSignature,</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">ExpectedPathSchema</FONT><BR ALIGN="LEFT"/>>];
    load [label=<<B>load_project()</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">checks each requirement;</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">select_signature picks a signature</FONT><BR ALIGN="LEFT"/>>];
    output_paths [label=<<B>get_output_paths()</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the work/ paths of</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the picked signature</FONT><BR ALIGN="LEFT"/>>];
    upstream_paths [label=<<B>get_upstream_paths()</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">the consumer's inputs, with the keys</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">of the field's ExpectedPathSchema</FONT><BR ALIGN="LEFT"/>>];
    upstream -> load [label="signatures"];
    field -> load [label="requirement"];
    load -> output_paths [label="picked signature"];
    output_paths -> upstream_paths;
}
```

### Signatures

A step lists what its tasks produce in `produces`.
If that depends on the task's config or its inputs, the step overrides the classmethod `task_produces(config, upstream)`.
For example, `ngs_mapping` produces `alignments [rna]` for STAR, and `variant_filtration` adds `filtered` to the tags of its input.
`variant_annotation` adds `annotated`, its tool and the genome build (`GRCh37`, `GRCh38`).
The steps after it keep these tags, and the calls report reads them for its build-dependent linkouts.

Each `depends_on` field states what it requires with a `DataSignature` in its `Annotated` metadata.
This annotation is the only place a requirement is declared.

A task can produce several signatures.
mutect2 calling produces all calls and the PASS calls, which carry the `filtered` tag.
`select_signature` picks the one a consumer reads: the `filtered` one when several satisfy the requirement.
So consumers get the PASS calls by default, and a requirement with `-filtered` gets all calls.

`get_output_paths` receives the picked signature and returns the matching file.
`get_upstream_paths(field, signature=...)` narrows the requirement for a single lookup.

### Checks and lookups

`load_project()` computes the signatures of every task in dependency order and checks each requirement before Snakemake builds any rule.
It also checks that the provider's `get_output_paths` for the picked signature has every required key of the field's `ExpectedPathSchema`.
So `get_output_paths` returns only the keys the task's config really produces.
A mismatch fails with a message that names both tasks.

Consumers read upstream files only through `get_upstream_paths(field, library_name=...)`.
It returns the provider's named outputs (`.bam`, `.vcf`, ...) wrapped in the field's `ExpectedPathSchema`.
The schema receives only the keys it declares.
Consumers never build a provider's paths themselves (`test_upstream_paths_come_from_contracts`).

`get_output_paths(config, signature, ...)` receives the config of the task whose outputs are requested.
It does not check the signature again.
It uses the signature only to choose between several outputs.

A field typed `list[str]` takes several tasks of the same kind, each checked like a single one.
The step lists them with `upstream_task_names(field)` and reads each with `get_upstream_paths(field, upstream_task=...)`.

### Granularity

Some steps read one result per pedigree, under the pedigree's name (`group_by: related`).
The `depends_on` field of such a step adds `RELATED_LEVEL` (`Granularity`) to its metadata.
A field that reads the result of a library's group at either level adds `ANY_LEVEL`.
Its step looks the result up with `upstream_entity()`.
`load_project()` checks the upstream task's resolved `group_by`.
For any other field, it rejects a per-library task behind one that writes per related group.
A lookup of a related group that the upstream task does not write fails when the DAG is built (`_check_upstream_library`).
The message names the fix.

### The `depends_on` vocabulary

`depends_on` keys are named after the data they bring in, not after the step that produces it.
The vocabulary has these keys: `reads`, `alignments`, `variants`, `copy_number`, `structural_variants`, `fusions`, `expression`, `hla_types`, `strandedness`, `qc`, `panel_of_normals`, `access`, `reference`, `features`, `proteome`, `dbsnp`, `common_variants`, `germline_resource`, `regions`, `excluded_regions`, `mappable_regions`, `par_regions`, `microsatellites`, `snp_loci`, `gene_regions`, `exons`, `kits`, `index`, `database`, `gcnv_model` and `umis`.
A step with two inputs of one kind qualifies them by role: `somatic_variants`, `germline_variants`, `combined_variants`, `annotated_variants`, `phased_variants` and `candidates`.
A step that needs only one input of a kind uses the plain key, such as `alignments`.
It narrows the input by the tags of the key's `DataSignature` instead (`alignments [rna]`).
`test_depends_on_keys_use_the_vocabulary` checks this.

### Required and optional keys

A key that every tool of a step reads is required and has no default.
A key that only some tools or settings read defaults to `""`.
It names its readers in its `Annotated` metadata (`snappy_pipeline/models`):

- `RequiredBy(...)` for those that need it,
- `ReadBy(...)` for those that read it if it is set.

A reader is a set of tools (none: any tool) and an optional condition on a setting:

- `when` is a dotted path into the step's config.
  The condition is on if every toggle along the path is enabled and every section on it is present.
  A bool at its end must be true as well.
- `value` compares the setting instead.
- `contains` tests membership in a list.

Examples:

- `RequiredBy(Tool.gatk4_hc_gvcf, ...)` for `dbsnp`,
- `RequiredBy(Tool.mutect2, when="mutect2.contamination")` for `common_variants`,
- `RequiredBy(when="intervals", value=Intervals.targeted)` for the kits of `gcnv_model`.

Several markers combine with "or".
For example, `kits` of `ngs_data_qc` has `RequiredBy(Tool.alfred)` and `RequiredBy(Tool.picard, when="picard.programs", contains=PicardProgram.CollectHsMetrics)`.
A validator of `SnappyStepModel` requires the key while a `RequiredBy` reader is active.
It rejects the key while no reader is active.
It reports all problems of a task in one `DependsOnError`.
Steps write no validators of their own for this.

Other optional keys default to `""`, meaning "not used".
No key defaults to a task name.

### Files and paths

A file that a step reads is a `depends_on` key if it must match the reference or other data of the project.
This holds even if only one tool reads the file.
Examples are regions, variants, a genome, a gene annotation and a mapper index (type `index`, tagged with the tool and the molecule).
Such a field has a `DataSignature` and an `ExpectedPathSchema`, and the step reads it with `get_upstream_paths()`.
Some files belong to the setup of one tool and do not depend on the reference, such as a jar or a script.
Such a file is a path in the tool's section, and the step passes it to the rule as an input.

A field that names a file or directory has the type `ConfigPath` (also in `list[...]` or `dict[str, ...]`).
`load_project()` expands its `pathvars` and resolves it against the project directory.
`snappy check-files` checks that the file or directory exists.

### Reference data and shared field types

A rule that reads the genome adds `self.parent.reference_inputs()` to its inputs.
These are the FASTA (`reference`) with its `.fai` (`reference_fai`) and `.dict` (`reference_dict`), named by the convention in `ExpectedReference`.
A rule that reads a VCF resource (`dbsnp`, `germline_resource`, `common_variants`) adds `self.parent.vcf_inputs(field)`.
These are the VCF (`<field>`) and its index (`<field>_index`), named by the convention in `ExpectedVcf`.
A BED file that GATK queries by interval, such as the mappable regions of gCNV, uses the schema `ExpectedIndexedRegions`.
The rule reads its `index` as well.

These inputs use the shared field types from `snappy_pipeline.workflows.abstract.protocol`:

| Input                                  | Field type         |
| -------------------------------------- | ------------------ |
| `reference` (genome FASTA)             | `Reference`        |
| `features` (gene annotation GTF)       | `Features`         |
| `dbsnp`                                | `Dbsnp`            |
| `germline_resource` (tag `population`) | `GermlineResource` |
| `common_variants` (tag `common`)       | `CommonVariants`   |
| the DNA alignments of a library        | `DnaAlignments`    |

The module also holds the shared schemas `ExpectedAlignments` (`bam`, `bai`) and `ExpectedVariants` (`vcf`, `vcf_tbi`).
Tags compose: one task with `tags: [population, common]` can serve both `germline_resource` and `common_variants`.
These keys usually name `external_data` tasks.

A step whose result is a bgzipped VCF declares its four output keys with `vcf_outputs(prefix)` from `snappy_pipeline.workflows.common`.
Its `get_output_paths` formats the same table as its parts.

### Databases and indexes

A tool's own database comes from `depends_on.database`: type `database`, tagged with the tool.
The field requires the tag `TOOL`, a placeholder that the loader replaces with the task's tool (`SnappyStepModel.tool_tag()`).
`ExpectedPathSchema.by_tool(...)` names the contract keys per tool.
`ngs_mapping` requires its index the same way.
Mapper and other indexes come from `depends_on.index` (type `index`, tags tool and molecule).

### Regions and kits

BED files of regions have the data type `regions` and the key `bed`.
The shared field type `Regions` accepts any of them, and `Exons` requires the tag `exons`.
A step-specific kind of regions is an `Annotated` field type in the step's model, with its own tag.
Add a field type with its own tag for a new kind of regions, so that a consumer cannot read the wrong kind.

`kits` uses the shared field type `Kits` (type `regions`, tag `kits`).
It names an `external_data` task with `kits`, which gives the `targets` and `baits` BED files of each enrichment kit.
A targeted library (`WES`, `Panel_seq`) uses the kit named exactly like its `library_kit`.
Other libraries have no kit (`BaseStep.library_kits`).
Long-read libraries, for example, use `libraryKit` for their chemistry.
`library_kit_paths(field, library_name)` returns the files of a library's kit.

A step calls `check_library_kits(field, library_names)` when it is set up.
So a kit missing from the kits task, often a typo, fails before Snakemake builds rules.
One message names all libraries with a missing kit.
Steps that need a kit for each library they process require one.
Steps whose per-library output is optional (the coverage report, a targeted gCNV model) skip libraries without a kit (`require=False`).
`common_kit()` returns the one kit of a group of libraries, for jobs over several libraries such as a panel of normals.

### `reads`

`reads` names an `external_data` task with the FASTQ files of the libraries, or a task that writes FASTQs (such as `adapter_trimming`).
An `external_data` task lists the files of each library in files tables (`library_files`; `FilesTable` in `snappy_pipeline/reads.py`), read once per run through `Project.files_table()`; `BaseStep.read_groups()` takes the read groups of a library from them.
`snappy scan-files` writes the tables with `ReadDiscovery` from the task's `scan` section (`search_paths`, `search_patterns`); nothing searches directories at run time.

### `FILTER` IDs

A step whose tasks write `FILTER` IDs other than their task name returns them from the classmethod `filter_ids(config)`.
Examples are the filters of mutect2 or the expression names of `vembrane` in tag mode.
`load_project()` rejects a task named `PASS`, and a task named like an ID of another task.
It also rejects a task that writes an ID that a task upstream of it writes.

### Configuration errors

A problem that the user has to fix in `config.yaml` or the sample sheets raises `ConfigError` (`snappy_pipeline.base`, a `ValueError`).
Its message states the problem only, such as `depends_on.reads names a task twice`.
The sheet loader names the sheet file.

The orchestrator attaches the task with `in_task(task)`.
It does so around its checks of one task and around `task_produces`.
It also does so around the constructors of the workflow classes (`check_config`, kits, library selection, relationships, sheet columns).
The error then reads `Task 'mapping' (step ngs_mapping): invalid configuration:` followed by one line per field.
The main Snakefile and `snappy check-files` print the error and exit.
Otherwise, Snakemake would frame it with its location and traceback lines (`test_run_prints_setup_errors_without_snakemake_framing`).

A `ValueError` or `TypeError` marks a programming error, such as a wildcard that names no library or a step part that returns no function.
Snakemake always frames an error raised in an input function with the rule and its wildcards.
So a check that can run at setup runs in the constructor.

## Paths

- Only truly variable values are wildcards: libraries, samples, tumor and normal, related groups (`related_group`: a pedigree or donor, named after it).
- The tool is part of the task configuration and never appears in a path.
  It is neither a directory prefix (`work/cnvkit.{library_name}/`) nor part of a file name.
- Files are named `{entity}.{suffix}`, with a dot after the wildcard, for example `{tumor_library}.gene_log2.txt`.
  Underscores stay inside multi-word suffixes.

`test_paths_follow_the_conventions` checks the last two rules on all DAG snapshots.
A value of the step's `tool` field counts as a tool name if it equals a run of whole path parts.
Path parts are separated by `/`, `.`, `_` and `-`, and the comparison drops these separators.
So `eb_filter` counts as the tool `ebfilter`.

## I/O footprint on HPC

On shared file systems such as cephfs, and with slurm, three counts cost more than the bytes moved: small files, metadata operations (stat, readlink) and jobs.
Keep all three low:

- no extra jobs for bookkeeping,
- no companion files without a consumer,
- no declared outputs nobody needs.

## Testing

`tests/snappy_pipeline/test_generated_configs_dryrun.py` builds the DAG of every generated task closure.
It compares each job (rule, wildcards, input, output and log paths, params, threads, resources and wrapper) with a snapshot in `tests/snappy_pipeline/snapshots/dag/`:

- A refactoring must leave the snapshots unchanged.
- An intended change is reviewed as a snapshot diff.
  Refresh the snapshots with `SNAPPY_UPDATE_SNAPSHOTS=1` (see {ref}`contributing`).
- Every closure must target at least one file.
  `EXPECTED_EMPTY` lists the exceptions with a reason.
  Closures that cannot be built yet are strict xfails in `KNOWN_BROKEN`.

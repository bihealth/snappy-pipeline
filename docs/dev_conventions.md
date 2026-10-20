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

## Step part API

The orchestrator loads the project once (`snappy_pipeline.orchestration.load_project`) and creates the workflow object of every task once.
A step Snakefile fetches it instead of building its own:

```python
from snappy_pipeline.orchestration import task_instance

wf = task_instance(config["__task_name__"])
```

`wf.get_task_config(name)` returns the task's own config (`name` is its step or task name).
It also returns the config of the upstream task set in `depends_on.<name>`.
Anything else is an error.

Snakefiles stay thin and call the workflow object for everything.
These are the methods a step part provides.
On the workflow object (`wf`), each method takes the step part name before the action, for example `wf.get_input_files("bwa", "run")`.

| Method                       | Returns                                                                    | Rule section                                                                            |
| ---------------------------- | -------------------------------------------------------------------------- | --------------------------------------------------------------------------------------- |
| `get_input_files(action)`    | a function of wildcards, returning a dict (named inputs), a list or a path | `input: unpack(wf.get_input_files(...))` for dicts, otherwise `wf.get_input_files(...)` |
| `get_params(action)`         | a function of wildcards (optionally also of `input`), returning a dict     | `params: args=wf.get_params(...)`                                                       |
| `get_output_files(action)`   | static paths (Snakemake does not accept functions here)                    | `output: **wf.get_output_files(...)`                                                    |
| `get_log_file(action)`       | static paths                                                               | `log: **wf.get_log_file(...)`                                                           |
| `get_resource(action, name)` | a function of wildcards, input, threads and attempt                        | `threads:` and `resources:`                                                             |

`BaseStep` raises a `TypeError` at parse time if `get_input_files` or `get_params` returns anything but a function.

A step part implements these through per-action methods that `BaseStepPart` dispatches to after validating the action.
These are `_get_input_files_<action>(wildcards)`, `_get_params_<action>(wildcards)`, `_get_output_files_<action>()` and `_get_log_file_<action>()`.
One `_get_log_file(action)` can replace the log methods of all actions.
Resource usage is declared, not computed.
`default_resource_usage` holds it for all actions and `resource_usage` per action, each as a property when it depends on the config.
`published_keys` names the output keys whose files are published, for all actions or by action.
`BaseStepPart.get_result_files` derives their `output/` paths for each output entity (`entity_wildcard`, other wildcards via `_result_wildcards`).
`BaseStep.get_result_files` collects them from all step parts.
Override these only where the derivation does not fit.

Wrappers read the params dict as `snakemake.params.args`.
Snakemake does not allow `unpack()` in `params:`, so every rule passes it as the single key `args`.

Wrappers never read `snakemake.config`.
Every file a wrapper reads, such as the reference FASTA, is a rule input, and every other value is a param.
Only inputs order the jobs that produce them first and trigger reruns when they change.
`test_wrappers_get_the_params_they_read` checks that the jobs in the DAG snapshots pass every `args` key their wrapper reads.
`test_files_are_inputs_not_params` flags params that hold file paths, against a baseline that only shrinks.

Rules that run a wrapper declare `benchmark: wf.get_benchmark_file(part, action)`.
This is the path of the rule's log below `benchmarks/` (`tasks/<task>/benchmarks/` in a project).
Snakemake writes the wall time, memory and I/O of each job there.
These files calibrate the resource requests.
Rules cached between workflows (`cache:`) cannot have a benchmark.
`test_wrapper_rules_write_a_benchmark` and `test_logs_and_benchmarks_are_unique_per_job` check this.

A rule's log holds only its log file (plus a copy of the script, if the wrapper writes one).
Snakemake's metadata in `.snakemake/` records the code, params, inputs and conda environment of every output.
The pinned environments (`environment.linux-64.pin.txt`) fix the package versions.
Logs are not published, and `snappy run` archives them per run.
Only published files, those linked into `output/`, have an md5 file, and every published file has one.
SODAR needs one per uploaded file.
A step declares `<file>.md5` next to `<file>` in the outputs, and the wrapper base classes write it.
Wrappers never read md5 files.
A step's `get_result_files` lists the published files below `output/`, each with its md5 file.
Bgzipped VCF and BED files also come with `.tbi` and `.tbi.md5`.
They mirror files below `work/`.
One local job per task, `<task>_publish` (`snappy_pipeline/publish.py`), links them and writes `output/.manifest.tsv`.
No other rule writes into `output/`.
`get_output_paths` returns `work/` paths, so tasks read each other's `work/` files and never `output/`.
`test_only_the_publish_job_writes_into_output`, `test_no_job_reads_from_output`, `test_logs_hold_only_logs`, `test_logs_are_not_published`, `test_md5_files_only_for_published_files`, `test_published_files_come_with_md5_and_index` and `test_wrappers_read_no_md5_or_conda_files` check this.

Every wrapper environment is pinned.
Snakemake builds it from the exact package URLs in `environment.linux-64.pin.txt` next to `environment.yaml`.
After changing an environment file, pin it again with `uvx snakedeploy pin-conda-envs --conda-frontend conda <file>` (this needs conda).
A symlinked `environment.yaml` gets a symlink to its target's pin file.
The `conda-autobump` GitHub workflow updates and re-pins the environments weekly, with one pull request per environment.
`test_wrapper_envs_are_pinned` checks that each environment has its pin file and that the pin file holds the exact versions of the environment file.

A step part lists the output keys of scatter fragments per action in `temp_outputs`.
It lists those of per-sample intermediates that only later rules of the same task read in `intermediate_outputs`.
`BaseStep.get_output_files` marks fragments `temp()` always and intermediates when the task sets `delete_intermediates`.
Published files, contract outputs, files other tasks read and project-wide resources (installs, prepared references, pedigree files) are neither.
Every temp file needs a declared reader, otherwise Snakemake deletes it right after its job.
A job that reads a temp file's index implicitly declares the index as input.
The DAG snapshots are built with `delete_intermediates: true`.
`test_intermediates_are_temp` and `test_temp_files_are_intermediates` check them.

Outputs for the Snakemake report are declared per action in a step part's `report_outputs`.
It maps the output key to a `ReportOutput` (subcategory, caption, `htmlindex` for a directory).
`BaseStep.get_output_files` flags them with `report()`, in the task's category and labelled by the wildcards of the path.
Captions are RST files in `report/` next to the step module.

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
- Override `get_input_files` or `get_params` only if it does more than dispatching, and then still return a function.
- `@dictify` and `@listify` keep the decorated function's signature (`functools.wraps`), which Snakemake inspects.
- Paths returned by input functions are relative (`work/...`, or upstream paths from `get_upstream_paths()`).
  Snakemake adds the task prefix to function results just as it does to static patterns.

The exception is an input function that needs Snakemake's `checkpoints` object.
The Snakefile calls it as `wf.get_input_files(...)(wildcards, checkpoints)` from a Snakefile-level function (`gcnv_model`).
Each task is a Snakemake module with its own `checkpoints` object, so a task cannot use the checkpoint of an upstream task.
It declares a local checkpoint instead, whose input is the upstream checkpoint output and whose output lists what the upstream checkpoint produced (`sv_calling_*` with `tool: gcnv`).

### Libraries, donors and pedigrees

Step parts read the sample sheets only through the library DataFrame of the workflow object, never through the biomedsheets objects.
`build_library_dataframe()` has one row per library, with its role, donor, cohort, parents, sex and the extra sample sheet columns; `include_background=True` adds the libraries of background data sets.
`library_row(name)` returns the row of one library, which `library_selection` picks if several data sets have a library of that name; read per-library facts (data set, folder, kit) through it, not by indexing the dataframe by name.
`output_entities` and `cohort_members` give the entities a task writes, with `library_selection` and `group_by` applied.
`primary_dna_libraries()`, `pedigrees()`, `pedigree_members()`, `pedigree_index()` and `children_with_both_parents()` give the germline view: one DNA library per donor, grouped by pedigree.

(dev_contracts)=

## Contracts

A task's outputs are described by `DataSignature` objects: a data type plus tags such as `dna`, `somatic` or `filtered`.

- A step lists what its tasks produce in `produces`.
  If that depends on the task's config or its inputs, the step overrides the classmethod `task_produces(config, upstream)`.
  For example, `ngs_mapping` produces `alignments [rna]` for STAR, and `variant_filtration` adds `filtered` to the tags of its input.
- `depends_on` keys are named after the data they bring in, not after the step that produces it: `reads`, `alignments`, `variants`, `copy_number`, `structural_variants`, `fusions`, `expression`, `hla_types`, `strandedness`, `panel_of_normals`, `access`, `reference`, `features`, `dbsnp`, `common_variants`, `germline_resource`, `regions`, `excluded_regions`, `mappable_regions`, `par_regions`, `microsatellites`, `snp_loci`, `gene_regions`, `exons`, `kits`, `index`, `database` and `gcnv_model`.
  A step with two inputs of one kind qualifies them by role, such as `somatic_variants` and `germline_variants`.
  `test_depends_on_keys_use_the_vocabulary` checks this.
- A file that a step reads is a `depends_on` key if it must match the reference or other data of the project, even if only one tool reads it.
  Examples are regions, variants, a genome, a gene annotation and a mapper index (type `index`, tagged with the tool and the molecule).
  Such a field has a `DataSignature` and an `ExpectedPathSchema`, and the step reads it with `get_upstream_paths()`.
  A file that belongs to the setup of one tool and does not depend on the reference, such as a jar or a script, is a path in the tool's section, which the step passes to the rule as an input.
- A field that names a file or directory has the type `ConfigPath` (also in `list[...]` or `dict[str, ...]`).
  `load_project()` expands its `pathvars` and resolves it against the project directory, and `snappy check-files` checks that it exists.
- A rule that reads the genome adds `self.parent.reference_inputs()` to its inputs: the FASTA (`reference`) with its `.fai` (`reference_fai`) and `.dict` (`reference_dict`), named by the convention in `ExpectedReference`.
- `reference` (genome FASTA), `features` (gene annotation GTF), `dbsnp`, `germline_resource` (tag `population`) and `common_variants` (tag `common`) use the shared field types `Reference`, `Features`, `Dbsnp`, `GermlineResource` and `CommonVariants` from `snappy_pipeline.workflows.abstract.protocol`.
  Tags compose: one task with `tags: [population, common]` can serve both.
  They usually name `external_data` tasks.
- A tool's own database comes from `depends_on.database`: type `database`, tagged with the tool.
  The field requires the tag `TOOL`, a placeholder that the loader replaces with the task's tool (`SnappyStepModel.tool_tag()`), and `ExpectedPathSchema.by_tool(...)` names the contract keys per tool.
  `ngs_mapping` requires its index the same way.
  Mapper and other indexes come from `depends_on.index` (type `index`, tags tool and molecule).
- BED files of regions have the data type `regions` and the key `bed`.
  The shared field type `Regions` accepts any of them, and `Exons` requires the tag `exons`.
  A step-specific kind of regions is an `Annotated` field type in the step's model, with its own tag.
  Add a field type with its own tag for a new kind of regions, so that a consumer cannot read the wrong kind.
- `kits` uses the shared field type `Kits` (type `regions`, tag `kits`): an `external_data` task with `kits`, which gives the `targets` and `baits` BED files of each enrichment kit.
  A targeted library (`WES`, `Panel_seq`) uses the kit named exactly like its `library_kit`; other libraries have no kit (`BaseStep.library_kits`), since long-read libraries use `libraryKit` for their chemistry.
  `library_kit_paths(field, library_name)` returns the files of a library's kit.
  A step calls `check_library_kits(field, library_names)` when it is set up, so that a kit missing from the kits task (often a typo) fails before Snakemake builds rules, with all such libraries in one message.
  Steps that need a kit for each library they process require one; steps whose per-library output is optional (the coverage report, a targeted gCNV model) skip libraries without a kit (`require=False`).
  `common_kit()` returns the one kit of a group of libraries, for jobs over several libraries such as a panel of normals.
- A key that every tool of a step reads is required and has no default.
  A key that only some tools read defaults to `""` and is listed per tool in the step's `TOOL_DEPENDENCIES`, which `validators.require_tool_dependencies` checks.
  Other optional keys default to `""`, meaning "not used".
  No key defaults to a task name.
- `reads` names an `external_data` task, a task that writes FASTQs (such as `adapter_trimming`), or the reserved value `data_sets`, which searches each library's FASTQs below the search paths of its data set.
  Data sets and `external_data` reads tasks share `LibraryFolders` (`search_paths`, `search_patterns`, `mixed_se_pe`) and one code path in `BaseStep.read_groups()`; only the folder name differs (the sheet's `folderName`, or the library name).
  No task may be called `data_sets`.
- Each `depends_on` field states what it requires with a `DataSignature` in its `Annotated` metadata.
  This annotation is the only place a requirement is declared.
- A task can produce several signatures. mutect2 calling produces all calls and the PASS calls, which carry the `filtered` tag.
  `select_signature` picks the one a consumer reads: the `filtered` one when several satisfy the requirement.
  PASS calls are therefore the default, and `-filtered` asks for all calls.
  `get_output_paths` receives the picked signature and returns the matching file.
  `get_upstream_paths(field, signature=...)` narrows the requirement for a single lookup.
- `load_project()` computes the signatures of every task in dependency order and checks each requirement before Snakemake builds any rule.
  It also checks that the provider's `get_output_paths` for the picked signature has every required key of the field's `ExpectedPathSchema`.
  So `get_output_paths` returns only the keys the task's config really produces.
  A mismatch fails with a message that names both tasks.
- Consumers read upstream files only through `get_upstream_paths(field, library_name=...)`.
  It returns the provider's named outputs (`.bam`, `.vcf`, ...) wrapped in the field's `ExpectedPathSchema`, and the schema receives only the keys it declares.
  Consumers never build a provider's paths themselves (`test_upstream_paths_come_from_contracts`).
- `get_output_paths(config, signature, ...)` receives the config of the task whose outputs are requested.
  It does not check the signature again.
  It only uses it to choose between several outputs.

## Paths

- Only truly variable values are wildcards: libraries, samples, tumor and normal, cohorts.
- The tool is part of the task configuration and never appears in a path, neither as a directory prefix (`work/cnvkit.{library_name}/`) nor in a file name.
- Files are named `{entity}.{suffix}`, with a dot after the wildcard, for example `{tumor_library}.gene_log2.txt`.
  Underscores stay inside multi-word suffixes.

`test_path_conventions_only_improve` checks the last two rules on all DAG snapshots against a baseline of existing violations (`tests/snappy_pipeline/snapshots/path_convention_violations.json`).
New violations fail the test, and fixed ones must be removed from the baseline.

## I/O footprint on HPC

On shared file systems such as cephfs, and with slurm, the number of small files, metadata operations (stat, readlink) and jobs costs more than bytes moved.
Keep all three low: no extra jobs for bookkeeping, no companion files without a consumer, no declared outputs nobody needs.

## Testing

`tests/snappy_pipeline/test_generated_configs_dryrun.py` builds the DAG of every generated task closure.
It compares each job (rule, wildcards, input, output and log paths, params, threads, resources and wrapper) with a snapshot in `tests/snappy_pipeline/snapshots/dag/`:

- A refactoring must leave the snapshots unchanged.
- An intended change is reviewed as a snapshot diff.
  Refresh the snapshots with `SNAPPY_UPDATE_SNAPSHOTS=1` (see {ref}`contributing`).
- Every closure must target at least one file.
  Exceptions are listed with a reason in `EXPECTED_EMPTY`.
  Closures that cannot be built yet are strict xfails in `KNOWN_BROKEN`.

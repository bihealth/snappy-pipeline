# Migrating from snappy 0.6

This page collects what changed between snappy 0.6 (the `main` branch) and the task-based snappy on Snakemake 9.
The pipeline configs in `.tests/test-workflow/pipelines/*/config.yaml` are complete examples of the new layout.

The page has four parts:

- The configuration model in detail, from {ref}`tasks <migrations_tasks>` to the {ref}`developer contracts <migrations_contracts>`.
- The new command, the changed defaults and the tooling, from {ref}`the snappy command <migrations_command>` to {ref}`CI <migrations_ci>`.
- {ref}`What was removed <migrations_removed>`, {ref}`what did not change, and what changes results on purpose <migrations_results>`.
- Step-by-step {ref}`migration notes <migrations_notes>` for existing projects, and an {ref}`example <migrations_example>` that translates a neoepitope config from main.

(migrations_tasks)=

## Task-based configuration

Instead of one block per workflow step, `config.yaml` lists **named tasks**.
Each task names a step type (`ngs_mapping`, `variant_calling`, ...) and carries its own config: the tool, the library selection, the dependencies and the tool settings.
Tasks form a DAG through their `depends_on` fields.
`external_data` tasks provide existing files, such as the genome, enrichment kit regions, indexes and databases.

By default `snappy run` builds every task and links its result files into `output/`.
`snappy run --task <name>` builds one task.
Adding `--frozen` makes it use the existing outputs of the upstream tasks instead of building them.

## Single tool per task

Workflows used to register all available tool implementations at once.
Now each workflow has a single `tool` field, an enum that selects exactly one implementation.
The tool's settings go into a section named after the tool:

```yaml
config:
  tool: gatk4_hc_gvcf
  gatk4_hc_gvcf:
    window_length: 10000000
```

This replaces main's `tools` lists.
It simplifies the Snakefiles, makes config validation more precise, and avoids ambiguity about which tool's resources or output paths apply.
Rarely, a step needs several tools that cannot be split into separate tasks.
The step's code then runs them, while the task's config still names one `tool`.

In the {ref}`example <migrations_example>` at the end of this page, the tool lists of main's `hla_typing` and `somatic_neoepitope_prediction` become five tasks:

```{graphviz}
:alt: The tool lists of main's hla_typing and somatic_neoepitope_prediction configs as five snake9 tasks

digraph tool_lists {
    rankdir=LR; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica"; labelloc=b;
    label=<<FONT POINT-SIZE="7" COLOR="#898781">Dashed arrows: a step config of main becomes tasks; solid arrows: depends_on.hla_types</FONT>>;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0];
    subgraph cluster_main {
        label=<<FONT POINT-SIZE="7" COLOR="#898781">main: step_config</FONT>>; color="#e1e0d9"; style="rounded"; labelloc=t;
        "main_hla" [label=<<B>hla_typing</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tools: {dna: [optitype], rna: [optitype, arcashla]}</FONT><BR ALIGN="LEFT"/>>];
        "main_neo" [label=<<B>somatic_neoepitope_prediction</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tools: [pvacseq, pvacfuse, pvacsplice]</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">tools_hla_typing</FONT><BR ALIGN="LEFT"/>>];
    }
    subgraph cluster_tasks {
        label=<<FONT POINT-SIZE="7" COLOR="#898781">snake9: tasks</FONT>>; color="#e1e0d9"; style="rounded"; labelloc=t;
        "hla_optitype_dna" [label=<<B>hla_optitype_dna</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">hla_typing, tool: optitype</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">DNA libraries</FONT><BR ALIGN="LEFT"/>>];
        "hla_optitype_rna" [label=<<B>hla_optitype_rna</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">hla_typing, tool: optitype</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">RNA libraries</FONT><BR ALIGN="LEFT"/>>];
        "hla_arcashla" [label=<<B>hla_arcashla</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">hla_typing, tool: arcashla</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">alignments of mapping_rna</FONT><BR ALIGN="LEFT"/>>];
        "neoepitopes" [label=<<B>neoepitopes</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">somatic_neoepitope_prediction, tool: pvactools</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">pvacseq, pvacfuse</FONT><BR ALIGN="LEFT"/>>];
        "neoepitopes_splice" [label=<<B>neoepitopes_splice</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">somatic_neoepitope_prediction, tool: pvactools</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">pvacsplice</FONT><BR ALIGN="LEFT"/>>];
    }
    "main_hla" -> "hla_optitype_dna" [style=dashed, color="#898781"];
    "main_hla" -> "hla_optitype_rna" [style=dashed, color="#898781"];
    "main_hla" -> "hla_arcashla" [style=dashed, color="#898781"];
    "main_neo" -> "neoepitopes" [style=dashed, color="#898781"];
    "main_neo" -> "neoepitopes_splice" [style=dashed, color="#898781"];
    "hla_optitype_dna" -> "neoepitopes" [color="#e87ba4b3"];
    "hla_optitype_rna" -> "neoepitopes" [color="#e87ba4b3"];
    "hla_arcashla" -> "neoepitopes" [color="#e87ba4b3"];
    "hla_optitype_dna" -> "neoepitopes_splice" [color="#e87ba4b3"];
    "hla_optitype_rna" -> "neoepitopes_splice" [color="#e87ba4b3"];
    "hla_arcashla" -> "neoepitopes_splice" [color="#e87ba4b3"];
}
```

## Library selection

Every task can filter which libraries it processes with a pandas `DataFrame.query()` expression:

```yaml
library_selection: "role == 'tumor' and extraction_type == 'dna'"
```

The expression runs against a DataFrame with one row per NGS library.

### Columns

The columns of the sample sheets, plus the derived ones:

| Column            | Description                                                                     |
| ----------------- | ------------------------------------------------------------------------------- |
| `library_name`    | The library, as the sheet names it                                              |
| `donor_name`      | The donor                                                                       |
| `sample_name`     | The sample, where the sheet has it                                              |
| `extraction_type` | `dna` or `rna`, from `library_type`                                             |
| `role`            | `tumor` for a library of a tumor sample (`is_tumor`), else `normal`             |
| `is_primary`      | The primary library among a donor's libraries of one role and extraction type   |
| `is_index`        | The index donor of its family                                                   |
| `sex`             | `male`, `female` or `unknown`                                                   |
| `family_id`       | The pedigree, where the sheet has it                                            |
| `disease_state`   | The PED phenotype of `is_affected`: 2 affected, 1 unaffected, 0 unknown         |
| `father_name`     | The father's primary normal DNA library, or `0`                                 |
| `mother_name`     | The mother's primary normal DNA library, or `0`                                 |
| `is_background`   | `true` for a library of a reference cohort                                      |

Every other sheet column is there under its own name, such as `library_type`, `library_kit`, `purity` or `batch_no`; see the sample sheets section of the overview.

Columns that are **optionally available** come from relationships.
Each `RelationshipDefinition` in the step's `default_relationships` or in the task's `relationships` adds a column.
The task's `relationships` override the step's defaults of the same name.

### Which libraries a task takes

A task's libraries come from three layers:

- What its step can process at all.
  This layer always applies.
  Examples are RNA libraries for expression quantification and fusion calling, tumor DNA libraries for mutect2, normal DNA libraries for the germline callers and for a panel of normals.
- Without `library_selection`, the step's default: the set its page names (mapping and `ngs_data_qc` take all libraries, the RNA steps all RNA libraries, the germline callers the normal DNA libraries, a panel of normals all normal DNA libraries), else the libraries of the task behind its first `depends_on` key that names a task, else `extraction_type == 'dna'`.
  The per-kind defaults of the data sets (tumor DNA on cancer sheets) are gone; a candidates task for somatic calls selects the tumors itself (`library_selection: "role == 'tumor'"`).
- `library_selection` replaces that default, for example to select a batch (`batch_no == 'run2'`) or to drop samples.
  It never goes beyond what the step can process.

### Background libraries

A library of a reference cohort has `is_background: true` in the sample sheet.
Such libraries are no task's results.
Mapping and panels of normals process them, and the SV and gCNV callers add them as cohort members or model samples.
`background_selection` narrows them per task with a query such as `batch_no == 'cohort_1000g'`.

## Group-by for aggregation level

Some steps (variant calling, annotation, filtration) can write one output per library or one per group of related libraries:

- `group_by: family`: one output per family (`family_id`), named after it.
  For a germline trio, this is one joint VCF for the whole family.
- `group_by: donor`: one output per donor (`donor_name`), which gives one joint VCF per donor (tumor and normal together).
- `group_by: related`: `family` when the task's libraries have a `family_id`, `donor` when none has.
- `group_by: library`: one output per library.
- Unset: the step's default, else the setting that upstream tasks declare, else `library`.
  Germline variant calling defaults to `related` and calls each pedigree jointly.

The setting propagates **downstream**: a filtration step after joint variant calling writes per pedigree too, unless it sets `group_by` itself.
So the aggregation level stays the same along a chain of dependent tasks.
When the project loads, steps that read one VCF per pedigree (VarFish export, peddy, IGV sessions, pedigree phasing) check that their upstream task has `group_by: related`.

(migrations_dependencies)=

## Dependencies and task-local files

A task gets its input files in one of two ways.
A **dependency** comes from another task, named in `depends_on`.
It has a data type and tags.
When the project loads, the consumer checks the data type, the tags and the keys of the files it reads.

Reference data are dependencies too: the genome, the gene annotation, dbSNP, enrichment kits, exons and other regions, and the mapper indexes.
`external_data` tasks provide existing files, and `reference_index` tasks build indexes.
So each reference file is named once, and a switch to another genome build only changes these tasks.

A file that only one tool reads can still be a dependency.
A bwa index has a format that only bwa reads, but it must match the genome that the other tasks use.
So it is `depends_on.index`, from a task with `produces: {type: index, tags: [bwa, dna]}`.

The type is the file contract (the keys a consumer reads), the tags say what the content is.
FASTQ read groups are `reads` (an `external_data` task with files tables, or `adapter_trimming`, which tags them `trimmed`), a FASTA is `sequences` (the genome: `{type: sequences, tags: [reference, dna]}`, a proteome: `[proteome]`), a GTF is `features`, and every BED file is `regions`, cnvkit's accessible regions included (`{type: regions, tags: [access]}` with the key `bed`).
There is no `raw` type.

A **task-local file** belongs to the setup of one tool and does not depend on the reference.
Examples are a jar or a Python file with vembrane helpers.
It is a path in the tool section, which the step passes to the rule as an input.

```yaml
# Entries of tasks:
- name: blacklist
  step: external_data
  config:
    produces: {type: regions, tags: [blacklist]}
    files: {bed: /path/to/blacklist.bed}

- name: blacklist_filter
  step: variant_filtration
  config:
    depends_on: {variants: calls, regions: blacklist}  # a dependency
    tool: regions
    regions: {action: exclude, mode: tag}

- name: lab_filter
  step: variant_filtration
  config:
    depends_on: {variants: blacklist_filter}
    tool: vembrane
    vembrane:
      mode: filter
      expression: 'is_reportable(QUAL, INFO["DP"])'
      context_files: [filters/lab_rules.py]  # a task-local file
```

All files of the maintained tools that depend on the reference are dependencies now.
The external export steps (`variant_export_external`, `wgs_cnv_export_external`, `wgs_sv_export_external`) are left for later.
They still run varfish-annotator with Jannovar `.ser` files from path fields.

### Regions and variant resources

BED files of regions have their own data type `regions` with the key `bed`.
The regions filter accepts any of them in `depends_on.regions`.
Steps that need a certain kind of regions require a tag (`exons`, `kits`, `mappable`, `microsatellites`, ...).

Population variants have the data type `variants` and a tag, as dbSNP does.
`depends_on.germline_resource` requires the tag `population`, and `depends_on.common_variants` the tag `common`.
On main, mutect2 read its resources from `mutect2.panel_of_normals`, `germline_resource` and `contamination.common_variants`.
They are now `depends_on` fields of `variant_calling` with the same names.
The panel of normals comes from a `panel_of_normals` task or an `external_data` task with `{type: models, tags: [pon, mutect2]}` and the key `vcf`.

### Gene annotation and tool databases

Every tool that reads a GTF reads `depends_on.features`.
RNA-SeQC needs a GTF prepared for it, so it requires the tag `rnaseqc` as well.

A tool's own database comes from `depends_on.database`: a task of type `database`, tagged with the tool.
Examples are the VEP cache and the mehari databases.
The field requires the tag of the reading task's own tool, which `load_project()` fills in.
So a VEP cache wired to a mehari task fails when the project loads.
Variant annotation and the VarFish export run mehari 0.51, which reads only transcript databases built by mehari 0.47 or newer (mehari-data-tx 0.14 or newer).
One mehari database task serves both.

### Sidecar files

Every job that reads the genome also has its `.fai` and `.dict` as inputs, so a missing sidecar fails when Snakemake builds the DAG.
Likewise, a job that reads dbSNP or population variants has the VCF's index as an input: `.tbi` for a `.vcf.gz`, `.idx` for a plain `.vcf`.

## Relative path resolution

Every path in `config.yaml` is relative to the directory of the config file, which is the project directory.
This holds for sample sheets, search paths, `external_data` files and files in tool sections.
`snappy run` passes the config as `--configfile`, so this also holds when Snakemake runs in another directory.
The top-level `pathvars` section names roots for paths, used as in `<refs>/genome.fa`.
It works like Snakemake's `pathvars` and has the same defaults.
snappy does not check whether the files exist when the project loads, so dry-runs work without them.
`snappy check-files` lists the missing ones.

## Work and output directories

Each task has its own directory `tasks/<task>/`.
`work/` holds every file the jobs write and is the contract between tasks.
A task reads the `work/` files of its upstream tasks, never their `output/`.
`output/` holds only the result files.
For each task with result files, one local `<task>_results` job links them from `work/` and writes `output/.manifest.tsv`.

```{graphviz}
:alt: The jobs of a task write into work/, which the results job and the downstream tasks read; the results job links the result files into output/

digraph work_output {
    rankdir=LR; nodesep=0.12; ranksep=0.35;
    bgcolor="transparent"; fontname="Helvetica";
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b", margin="0.08,0.03"];
    edge [arrowsize=0.5, penwidth=1.0, color="#898781"];
    "jobs" [label=<<B>jobs of a task</B>>];
    "work" [style="rounded,dashed", label=<<B>tasks/&lt;task&gt;/work/</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">every file the jobs write</FONT><BR ALIGN="LEFT"/>>];
    "results" [label=<<B>&lt;task&gt;_results</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">one local job</FONT><BR ALIGN="LEFT"/>>];
    "output" [style="rounded,dashed", label=<<B>tasks/&lt;task&gt;/output/</B><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">links to the result files</FONT><BR ALIGN="LEFT"/><FONT POINT-SIZE="7" COLOR="#52514e">.manifest.tsv</FONT><BR ALIGN="LEFT"/>>];
    "downstream" [label=<<B>jobs of a downstream task</B>>];
    "jobs" -> "work";
    "work" -> "results";
    "results" -> "output";
    "work" -> "downstream";
}
```

Only result files get an `.md5` file.
A job's log is only its log file, without main's conda listings and `.md5` files of logs.
`snappy run` archives the logs of each run.
Scatter fragments are always temporary.
`delete_intermediates: [intermediates]` makes the per-sample intermediates temporary, too.
These are the files that only later rules of the task read.
`delete_intermediates: [bam]` does so for the BAM of a mapping task that a later task processes further.

(migrations_contracts)=

## Developer contracts: produces and DataSignature

For developers writing new workflow steps, `load_project()` validates the dependency chains before Snakemake builds any rule.

- **`produces`** declares the `DataSignature` values of a task's outputs, such as `ALIGNMENTS` with tag `dna`.
  Where the outputs depend on the task's config or its inputs, the classmethod `task_produces(config, upstream)` declares them instead.
- **`depends_on` fields** carry the `DataSignature` they require and an `ExpectedPathSchema` of the files they read.
- **`DataSignature`**: a frozen dataclass of a `DataType` enum and a frozenset of tags.
  A required tag is a string that must be present, or a tuple of which one must be present.
  A string with a `-` prefix names a tag that must be absent, such as `"-somatic"`.

`load_project()` checks that each upstream task produces what the `depends_on` field requires, and that an `external_data` task has the keys the consumer reads.
This catches mismatches such as wiring an RNA alignment step into a variant caller that expects DNA, or a germline caller into a somatic pipeline.
It also checks that an index from a `reference_index` task is built from the genome the reader names, and that `FILTER` IDs cannot mean two things.
Consumers read upstream files only through `get_upstream_paths()`.
On main, consumers built the upstream paths themselves and assumed that they would match.

(migrations_command)=

## The `snappy` command

A single `snappy` command with subcommands replaces the old entry-point scripts.
{ref}`usage_commands` lists all subcommands.

| main                   | snake9                                                               |
| ---------------------- | -------------------------------------------------------------------- |
| `snappy-start-project` | `snappy init`                                                        |
| `snappy-start-step`    | `snappy task add`                                                    |
| `snappy-snake`         | `snappy run`                                                         |
| `snappy-pull-sheet`    | `snappy pull-sheet`                                                  |
| `snappy-refresh-step`  | none: the executor comes from a Snakemake profile                    |
| `snappy-slurm-status`  | none: the SLURM executor plugin of Snakemake keeps track of the jobs |

`snappy pull-sheet` writes the batch of SODAR verbatim into the `batch_no` column, a string.
Main declared it an `integer` and wrote the first number of the batch, or 1 without a batch column.

### `snappy run` options

- `--task <name>`: build only the named task and the upstream files it needs, and link its result files (default: every task).
- `--frozen`: with `--task`, use the existing outputs of the upstream tasks and never rebuild them.
- `--directory`: the project directory (default: the current directory).
- `--verbose`: log the Snakemake command line and other debug messages.
- Everything after `--` is passed verbatim to Snakemake.

## Changed defaults

### Default workflow profile

snappy always applies its workflow profile (`snappy_pipeline/apps/profile/`).
It sets:

- `software-deployment-method: [conda, container]` with `sdm-container-runtime: apptainer`: jobs run in their conda environments, and the tools with a container image run in it with apptainer, without `--use-conda` or `--use-apptainer`.
- `logger: snkmt`: structured logging to the SQLite database `.snakemake/log/snkmt.sqlite`.
  `snappy watch` shows it live.
- `persistence-backend: db`: a single SQLite file for Snakemake's metadata instead of many small files in `.snakemake/`.
  This matters on shared network file systems such as CephFS and NFS, which handle one large file much better than thousands of small ones.

snappy passes this profile with `--workflow-profile`, as Snakemake 9 recommends.
The executor and its settings come from a Snakemake profile of the system or the user, selected with `--profile` after `--` or with `SNAKEMAKE_PROFILE`.

### Cluster execution

The executor comes from a Snakemake profile of the system or the user, not from snappy.
Pass it as `--profile <name>` after `--`, or set `SNAKEMAKE_PROFILE`.
The old `--profile-snappy-pipeline` flag, the `pipeline_job.sh` job script and the `snappy refresh` and `snappy status` commands are gone.
The Snakemake process runs in a terminal multiplexer on the login node or in an allocation (`salloc`), as the documentation of the SLURM executor recommends.
{ref}`usage` has an example profile for SLURM.

### Output path patterns

Result file paths are now much shorter.
Main put the tool name and the step's substage into every path, even with only one tool and one substage.
Examples are `bwa.filtered.vcf` and `bwa.raw.vcf`.
The task config already defines both, so the paths no longer repeat them.
Paths now contain only the parts that vary, usually the library name.
For related libraries, the path holds the name of the pedigree or donor instead.
Multi-tool workflows and substeps whose paths would otherwise collide still use more wildcards.
The default is a single wildcard.

`report/` holds what the Snakemake report shows.
Data for other tools is in `out/`, such as the coverage VCF and bigWig files of `bam_collect_doc` (main: `report/cov/`).
Tasks whose tools write statistics add one job that writes the result file `qc_summary/out/qc_summary.tsv`, with one row per library.
These are `ngs_mapping`, `ngs_data_qc` with `alfred`, `picard` or `ngs_chew`, `variant_calling` with `bcftools_stats`, and `variant_checking`.
The job also renders tables and charts of all libraries for the report with datavzrd.

Files that main's wrappers wrote without declaring them are declared outputs now:

- PureCN's plots: `{library}.pdf`, `.local_optima.pdf`, `.segmentation.pdf` and `.chromosomes.pdf` instead of `_local_optima.pdf` and so on.
- STAR's `{library}.Log.final.out`.
- Salmon's `{library}.lib_format_counts.json` and auxiliary directory `{library}.aux` (main: `out/aux/`).
  Salmon's own log is part of the job log.

### Rule names and the temporary directory

The rules of the wrapper actions are generated from the step parts and named `<task>_<step>_<part>_<action>`, for example `calling_variant_calling_mutect2_run`.
A rule that main named differently (`variant_annotation_vep`, `somatic_msi_calling_mantis_msi2`) now carries its action.
So `--until`, `--omit-from` and `--allowed-rules` arguments need the new names.

The temporary directory of the jobs is `resources.tmpdir` of the project config (main: `global_config.tmpdir`, with `$TODAY`), else `$TMPDIR`.

### MSI calling and popdel

`somatic_msi_calling` leaves out tumor libraries without a matched normal and warns about them, because MANTIS2 needs the normal.

`sv_calling_wgs` with `tool: popdel` calls with Snakemake's scatter-gather.
Main ran one job per window and read the windows from the `.fai` when Snakemake built the DAG.
The windows stay the same.
Each contig that `ignore_contigs` keeps is cut into windows of `popdel.window_size` bases (10 Mbp) from its start.
Each window is padded by `popdel.max_sv_size` bases (20 kbp) on both sides, within its contig.

The new `popdel.num_jobs` (24) is the number of jobs that share the windows in genome order.
Each job gets at least one window and about equal bases, and calls its windows one after the other.
A reference with fewer windows than jobs is an error.
Except at the start of a contig, popdel's region of a window starts one base later than on main, at the first base of the padded window.
Main passed the 0-based start as the 1-based one, so calls within `popdel.max_sv_size` of a window edge can differ from main.

## Snakemake 7 to Snakemake 9

The pipeline now requires Snakemake 9 or newer.
The most visible change is that Snakemake 9's **module** system replaces the **subworkflows**.
Each step's Snakefile is imported as a module with a `replace_prefix` mapping.
So each task has its own output and work directories, such as `tasks/sample_bwa/output/` instead of a flat `output/`.
This lets one project run the same step type several times without file collisions.

More wrappers use the typed wrapper classes `ShellWrapper` and `PythonWrapper`.
R code runs in heredocs of a `ShellWrapper`.
These improve type checking and make wrapper signatures more explicit compared to the plain `Wrapper` base class.

The move to Snakemake 9 also required updates of several dependencies to stay compatible with Snakemake's own dependencies.
Among them are pandas and pydantic.

snappy uses the development version of Snakemake, which deploys software through plugins.
The conda plugin builds the wrapper environments from their pin files, below `.snakemake/software/`.
`snappy run -- --deploy-software-envs` creates them without running jobs (formerly `--conda-create-envs-only`).

Some R packages are not on conda-forge or bioconda in the needed versions (for sequenza and scarHRD).
The post-deploy script of the wrapper's environment (`environment.post-deploy.sh`) installs them at fixed versions or commits when the environment is created.
PureCN and pVACtools run from their container images (`software: container(...)`), which Snakemake pulls and runs with apptainer.

## Project management with pixi

We adopted [pixi](https://pixi.sh) as the project manager.
<https://pixi.sh/latest/#installation> describes how to install it.
Snappy has conda dependencies for the system tools it calls itself, such as samtools, bcftools and htslib, and PyPI dependencies for the Python packages; the tools of the steps come from the wrapper environments.
Pixi handles both.

The pixi configuration defines several development tasks:

- `pixi run -e dev test`: run the test suite
- `pixi run -e dev lint`: run all linters (ruff check, ruff format, snakefmt, shfmt and the module layout, all in check mode).
  `mypy` is not in use.
- `pixi run -e dev fmt`: auto-format with ruff
- `pixi run -e dev snakefmt`: auto-format Snakemake files
- `pixi run -e dev srcfmt`: run all formatters and the module layout (`fmt`, `snakefmt`, `shfmt` and `layout`)

For developers, a single `pixi install` replaces the old `conda env create` and `pip install` steps.

## Build system and packaging

The build backend is `hatchling` instead of `setuptools`.
The old `setup.py`, `setup.cfg`, `environment.yml` and `Makefile` are removed.
Pixi manages the environments, and the build configuration is entirely in `pyproject.toml`.

(migrations_ci)=

## CI and test suite modernization

CI sets up its environments with `prefix-dev/setup-pixi` instead of the legacy `conda-cache` GitHub action.
The end-to-end tests (`.github/workflows/ci-e2e.yml`) run the test pipelines in `.tests/test-workflow/pipelines/` with `pixi run snappy run`.
Their test data is checked in with Git LFS.
They dry-run and then run the cancer, germline and reference pipelines on small test data.
Wrapper environments are pinned (`environment.linux-64.pin.txt`).
In a pull request, a changed environment file must come with its changed pin file (`pixi run -e pins update-pins <file>`).

(migrations_removed)=

## What was removed

### Workflows and tools

| Removed                                                                                                                                          | Reason or replacement                                                                                                                                         |
| ------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **`somatic_wgs_sv_calling` workflow**: a specialised SV calling pipeline for whole-genome data with its own complex config                       | It was not widely used. SV calling for targeted data with `sv_calling_targeted` is unaffected.                                                                |
| **Jannovar statistics report** (`jannovar_stats` of `variant_calling`)                                                                           |                                                                                                                                                               |
| **`mbcs` tool of `ngs_mapping`** and its nested Snakemake run (`nested_snakemake.py`), with `depends_on.common_variants` of `ngs_mapping`        | Molecular barcodes (UMIs) take an `umi_mapping` task with fgumi, BQSR an `alignment_processing` task (migration note 20).                                     |
| **`fusioncatcher` tool of `somatic_gene_fusion_calling`**                                                                                        | It had no wrapper, not even on main, so it could not run.                                                                                                     |
| **`defuse`, `hera`, `jaffa` and `pizzly` tools of `somatic_gene_fusion_calling`**                                                                | They were untested and probably broken. `arriba` (recommended) and `star_fusion` remain.                                                                      |
| **main's germline filters** (`variant_denovo_filtration`, the germline filters of main's `variant_filtration` and the `protected` rescue filter) | They no longer ran on main (a deleted shell library, Jannovar inputs). `vembrane` expressions in a `variant_filtration` task cover their per-record criteria. |

### Commands, profiles and jobs

| Removed                                                                                                                                                                                                                                                                  | Reason or replacement                                                                                                                     |
| ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------- |
| **Old `tpls/profile/` directory**                                                                                                                                                                                                                                        | Superseded by the workflow profile `snappy_pipeline/apps/profile/` and a Snakemake profile of the system or the user for the executor.    |
| **`--profile-snappy-pipeline` flag, `pipeline_job.sh`, `snappy refresh` and `snappy status`**                                                                                                                                                                            | Cluster execution is configured in a Snakemake profile.                                                                                   |
| **Legacy `snappy-snake` entry point**                                                                                                                                                                                                                                    | Use the `snappy` command instead.                                                                                                         |
| **The additional command line tools** `snappy-bed_filter_jaccard`, `snappy-fix_vcf`, `snappy-genome_windows`, `snappy-quickvenn`, `snappy-vcf_first_header`, `snappy-vcf_filter_denovo` and `snappy-ped_to_vcf_header`, with `add_ped_header.sh`                         | Nothing in the pipeline ran them. `snappy-vcf_filter_to_info` and `snappy-vcf_filter_from_info` stay for the read-backed phasing wrapper. |
| **`config.yaml.bak` copies** of `snappy task add` and `snappy init --task`                                                                                                                                                                                               | The commands print the change as a diff, and git keeps the history.                                                                       |
| **R package install jobs of sequenza and scarHRD** and their `work/R_packages/` directories                                                                                                                                                                              | The environments install these packages when they are created.                                                                            |
| **Container builds of PureCN and pVACtools**: the rules `panel_of_normals_purecn_install` (the `singularity` wrapper, `work/containers/`) and `somatic_neoepitope_prediction_pvactools_install` (`work/container/`), and the key `container` of a PureCN panel's outputs | Snakemake pulls the images.                                                                                                               |

### Config fields

**Multi-tool list configs**: the `tools:` list in step configs is no longer supported.
Use `tool:` with a single value.

These config fields are removed because nothing read them or they duplicated the upstream config.
A config that sets one of them fails validation:

- `cbioportal_export`: `vcf2maf.annotation_tool` (the export requires VEP annotations), `copy_number_alteration.copy_number_tool` and `patient_info`.
- `somatic_cnv_checking`: `cnv_assay_type`.
- `hla_typing`: `arcashla.mapper`.
- `somatic_targeted_seq_cnv_calling`: `somatic_variant_caller` of PureCN.

These fields are gone as well:

| Step or section                                        | Removed fields                                                                                                                                                                                         | Why                                                                                                     |
| ------------------------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ------------------------------------------------------------------------------------------------------- |
| The mutect2 sections                                   | `use_profile`, `restart_times`, `max_jobs_per_second`, `max_status_checks_per_second`, `debug_trunc_tokens`, `keep_tmpdir`, `job_mult_memory`, `job_mult_time`, `merge_mult_memory`, `merge_mult_time` | The settings of the removed nested Snakemake runs, which nothing read. `num_jobs` and `num_cores` stay. |
| `homologous_recombination_deficiency`                  | `scarHRD.length`, `depends_on.reference`                                                                                                                                                               | scarHRD read neither, only the `seqz` file of the copy number task.                                     |
| `somatic_targeted_seq_cnv_calling` with `tool: purecn` | `purecn.genome_name`, `purecn.seed`                                                                                                                                                                    | Its PureCN jobs take the values of the panel of normals task, as they take its image.                   |
| `somatic_targeted_seq_cnv_calling`                     | `purecn.enrichment_kit_name`                                                                                                                                                                           | No wrapper read it.                                                                                     |
| The `cnvkit` section of `panel_of_normals`             | `genemetrics`, `segmetrics`, `breaks`, `drop_low_coverage`                                                                                                                                             | The panel report never read them.                                                                       |

(migrations_results)=

## What did not change, and deliberate result changes

The rules run the same tools, such as BWA, STAR and GATK.
The per-step Snakefiles (`workflows/*/Snakefile`) still use the same patterns: `wf.config`, `wf.get_input_files()`, `wf.get_params()` and wrappers.

### Comparison with main 0.6.3

A real-run comparison of the e2e pipelines against main 0.6.3 gives the same results, with two exceptions:

- Main's bcftools 1.15.1 cannot parse the GT of the dkfz-filtered calls and drops every record.
  snake9's bcftools 1.23.1 keeps and tags them.
- `FILTER` IDs are task names (migration note 10).

### Deliberate result changes

Outside the e2e pipelines, results change on purpose where main was wrong or ignored a setting:

| Step or tool         | Change                                                                                                                 | Main                                                                                                           |
| -------------------- | ---------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------- |
| MELT                 | GroupAnalysis reads the discovery output of all pedigree members.                                                      | Read only the first member's.                                                                                  |
| mutect2 and sequenza | The `ignore_contigs` defaults apply (migration note 14).                                                               |                                                                                                                |
| arcasHLA             | Gets the configured `zygosity_threshold` (migration note 17).                                                          |                                                                                                                |
| VEP                  | Annotates multiallelic records split into one record per allele (`bcftools norm`).                                     |                                                                                                                |
| ebfilter             | The panel draws from all normal DNA libraries.                                                                         | Drew from the normals matched to a tumor.                                                                      |
| Control-FREEC        | Takes each donor's sex from the sample sheet.                                                                          | Used `XX` for every sample.                                                                                    |
| PureCN               | A panel reads the GenomicsDB of its `mutect2` panel.                                                                   | Failed there, because the `genomicsdb` R package of the PureCN image needs integers where GATK writes strings. |
| BWA mappers          | Read a fixed 100 Mbp per batch (`bases_per_batch`, `-K`), so their alignments no longer depend on `num_threads_align`. |                                                                                                                |
| BWA mappers          | Their read group IDs are the names of the discovered read groups.                                                      | `<library>.0`, `<library>.1`, ... for `bwa`, and one read group for all reads for `bwa_mem2`.                  |
| `bwa_mem2`           | Maps single-end reads.                                                                                                 | Read the `right` file, which single-end reads do not have.                                                     |
| `bwa_mem2`           | Fails when it cannot decompress the reads.                                                                             | Wrote an empty BAM.                                                                                            |
| BWA mappers and STAR | Mask the duplicates of single-end reads.                                                                               | samblaster aborted on reads without mates, and `bwa` could not read single-end files at all.                   |

### Pedigrees and joint germline calls

The family column of the PED files is the pedigree's name, the sheet's `family_id`.
The trio PED of a child who is not the index names the child's family.
The index of a family has a normal DNA library, and where several donors have one, the sheet marks the index with `is_index`; a family whose index had none fails the sheet check (main dropped its members silently).

A joint germline call takes one DNA library per pedigree member, so a donor's RNA libraries and further DNA libraries stay out of it.
On a cancer sheet, it takes the normal, or the primary tumor of a donor without a normal.
Before, it took every library of the donor, the tumor included.

(migrations_notes)=

## Migration notes

To migrate an existing project, go through these notes.

### Configuration layout and commands

1. Wrap each step's configuration in a `tasks:` list entry with `name` and `step` fields.

2. Replace `tools: [bwa]` with `tool: bwa`.
   Move tool-specific config out of the list into a key of the task config named after the tool.

3. Remove `--profile-snappy-pipeline` from scripts.
   snappy applies its workflow profile automatically.
   Put the cluster settings into a Snakemake profile of the system or the user.
   Pass it as `--profile <name>` after `--`, or set `SNAKEMAKE_PROFILE`.

4. Use `snappy run` instead of the old entry points.

5. snappy rejects the old `config.yaml` structure without a `tasks:` key.
   Run `snappy init` in a fresh directory to see the new layout.

### Reference files and sample sheets

6. Move files that depend on the reference out of tool sections into `external_data` tasks named in `depends_on` (see {ref}`migrations_dependencies`).
   For example, a VEP task now needs a cache task in `depends_on.database`.

7. snappy no longer looks for sample sheets in the parent directory of the project.
   `sample_sheets` lists them (default: `samplesheet.tsv` next to `config.yaml`); the paths of `config.yaml` are relative to the project directory, absolute, or below a `pathvars` root.
   The biomedsheets format is gone: a sheet is a plain TSV with one header line and one row per library, which `snappy sheet check` checks (see the sample sheets section of the overview).
   Convert an old sheet by hand or with a script: `patientName` is `donor_name`; `sampleName` is `sample_name`, as the full name (`P001-N1`); `isTumor` is `is_tumor` (`true`/`false`); `fatherName` and `motherName` are `father` and `mother`, the donor names, empty instead of `0`; `sex` is `male`/`female`/`unknown`; `isAffected` is `is_affected` (`true`/`false`); `familyId` is `family_id`; `libraryType` is `library_type`; `folderName` is `folder_name`; `libraryKit` is `library_kit`; `hpoTerms` is `hpo_terms`; the other custom fields are columns of the same name in snake case; `extractionType` is dropped; the library names are written out as biomedsheets generated them (`P001-N1-DNA1-WES1`), so the files keep their names.
   Mark the index donor of each family with `is_index` (biomedsheets took the first affected sequenced donor by name) and, where a donor has several libraries of one role and extraction type, the primary one with `is_primary`.
   The data set `type`, `naming_scheme`, `pedigree_field` and `is_background` are gone; `is_background` is a sheet column.

8. snappy searches no directories at run time.
   An `external_data` task lists the FASTQ files of each library in a files table (`library_files: [reads.tsv]`, the columns `library`, `readgroup`, `left`, `right`), and the tasks name it in `depends_on.reads` instead of `data_sets`.
   `snappy scan-files` writes the table from the files below `raw` in a folder named like the sheet's `folderName` (else the library) or named after it, with the common namings (`_R1_001.fastq.gz`, `.R1.fastq.gz`, `_R1.fastq.gz`, `_1.fastq.gz`, also `.fq.gz`); `--search-path` and `--pattern key=regex` (regular expressions for the path below the folder or the file name, whose `(?P<readgroup>...)` group pairs the mates) change that, and the table records the command.
   Run it once per project and again after adding libraries; `snappy run --scan` does it first.
   Existing alignments or variants enter the same way, as a table with the columns `library`, `bam`, `bai` or `library`, `vcf`, `vcf_tbi`.

### Steps and tools

9. VEP annotates all transcripts in a task of its own with `vep.pick: false`, instead of the `.full.vcf.gz` side output.
   Consumers that need all transcripts (`use_all_transcripts` of neoepitope prediction, `use_all_transcripts_for_*` of combine_variants) read such a task.
   VEP plugins (`vep.plugins`) are local `.pm` files.
   snappy no longer downloads them from a URL as main did.
   Annotations are tagged with the annotator (`vep` or `mehari`), so a consumer can require one.

10. The `bcftools`, `regions` and `ebfilter` filters of `variant_filtration` use the task name as FILTER ID.
    Main used `<filter>_<position in filter_list>`, such as `bcftools_2`.
    Expressions that test FILTER values need the new names.
    No task may be named `PASS` or like a FILTER ID that another task writes.
    These IDs include the filters of mutect2, `LowQual` of the GATK callers, the expression names of `vembrane` in tag mode, and `bPcr` and `bSeq` of `dkfz`.
    No task may write an ID that a task upstream of it writes, either.

11. Every tool's config section takes `extra_arguments`, a list of additional arguments for the tool's main command line.
    `variant_checking`, `repeat_expansion` and other single-tool steps without a tool section take it on their step model.
    `extra_arguments` replaces the former `extra_args` fields of `bwa`, `bwa_mem2`, pVACtools, `reference_index` and `vembrane` (whose value was a string).
    A tool with several commands passes them to one, such as HaplotypeCaller for `gatk4_hc_gvcf`, `segment` for cnvkit or `GermlineCNVCaller` for gCNV.
    The tumor phasing tools (`gatk_read_backed_phasing_tumor`, `whatshap_tumor`) and `somatic_variant_signatures` (`tabulate_vcf`) have new sections for it.

12. The BED file of mappable regions for gCNV needs its index next to it: `.idx` from `gatk IndexFeatureFile`, or `.tbi` for a bgzipped BED.
    `snappy check-files` lists the index if it is missing.

13. `somatic_neoepitope_prediction` is one task with `tool: pvactools`.
    `pvactools.pvacseq` (on by default), `pvacfuse` and `pvacsplice` are toggles instead of the `tools` list.
    Its inputs come from `depends_on` instead of `path_*` and `tool_*` fields:

    - `hla_types` lists the HLA typing tasks instead of `tools_hla_typing`.
    - `alignments` (RNA) replaces the mapper of `pileup`.
    - `pileup`, `quantification` (now `expression`) and `phasing` are pvacseq options.
    - `path_netchop` is `net_chop.path`.

    Phasing runs in a `variant_phasing` task, which the neoepitope task reads through `depends_on.phased_variants`.
    Its tool is `gatk_read_backed_phasing_tumor` or `whatshap_tumor`.
    `gatk_read_backed_phasing_tumor` runs GATK 3 ReadBackedPhasing on the combined variants with the tumor DNA, as on main.
    `whatshap_tumor` runs WhatsHap and writes ReadBackedPhasing's `HP` tag, which pVACseq reads.
    The phased variants should be annotated for all transcripts, because pVACseq matches nearby variants by the transcript of the somatic variant.

    The variants must come from VEP with the pVACtools plugins `Frameshift` and `Wildtype`.
    The task checks this when Snakemake builds the workflow.
    For phasing, the germline calls must be VEP-annotated, too.
    Predictions that read different annotations (pvacsplice with all transcripts, pvacseq without) are separate tasks.
    The proteome is an external FASTA (`depends_on.proteome`).
    Building it from germline variants is not ported.
    The {ref}`example at the end of this page <migrations_example>` translates a neoepitope config from main.

14. `ignore_chroms` is renamed to `ignore_contigs` in every step (`variant_calling`, `variant_phasing`, `sv_calling_wgs`, `panel_of_normals` and the `sequenza` section).
    The entries are regular expressions matched with `re.search`, where main's somatic steps used globs.
    `HLA-*` becomes `^HLA-`, and `*_decoy` becomes `_decoy$`.
    In `variant_calling`, the list is a field of each caller's section instead of the step.
    A YAML anchor can share one list between tasks.
    By default, mutect2 skips decoys, HLA, EBV and the hs37d5 extras, like the germline callers.
    On main, it called on every contig.

15. `bbduk.adapter_sequences` takes only reference files, with paths relative to the config file.
    BBDuk's keywords (`adapters`, `artifacts`, `phix`, `lambda`, `pjet`, `mtst`, `kapa`) go in `bbduk.builtin_references`.

16. The `tag_mode` of main's vembrane filter is gone, because an expression names the variants to tag.
    `tag_mode: include` becomes the negated expression, `not (...)`.

17. `arcashla.zygocity_threshold` is renamed to `arcashla.zygosity_threshold`, and it now reaches arcasHLA.
    Main never passed it, so arcasHLA ran with its own default 0.1 instead of the configured 0.15.

18. `somatic_wgs_cnv_calling` with `tool: cnvkit` needs `depends_on.panel_of_normals`, a whole-genome `panel_of_normals` task with `tool: cnvkit`.
    Without normals, that task builds a flat reference.

19. `purecn.container` and `pvactools.container` name the images that Snakemake pulls.
    By default, they are pinned by digest: `markusriester/purecn` of 2024-05-28 and `griffithlab/pvactools:7.1.3`.
    The path of a local `.sif` file works too.
    `pvactools.container: null` is no longer valid.
    A PureCN calling task runs the image of its panel of normals task.
    Jobs in a container see the project, `$HOME`, `/tmp` and the snappy installation.
    Add other directories that they read to `APPTAINER_BIND`, as {ref}`usage` describes.

### `ngs_mapping` features as tasks

20. `tool: mbcs` of `ngs_mapping` is removed with its AGeNT steps.
    Libraries with molecular barcodes (UMIs) take an `umi_mapping` task with `tool: fgumi` ({ref}`step_umi_mapping`), whose consensus BAM files the callers read like those of a mapping task; the read structure comes from the sheet's custom field `readStructure`.
    Its BQSR (`recalibrate`, `depends_on.common_variants`) is an `alignment_processing` task with `tool: gatk4_bqsr`.
    It follows a mapping task with the tool of `mbcs.mapping_tool`, as the recipe in {ref}`somatic_ngs` shows.
    With main's mbcs output without BQSR as input, `gatk4_bqsr` writes main's recalibrated BAM.

21. `ngs_chew_fingerprint` of `ngs_mapping` becomes an `ngs_data_qc` task with `tool: ngs_chew` and `depends_on.alignments` naming the mapping task.
    Its `extra_arguments` go to `ngs_chew.extra_arguments`.
    The fingerprints are `{library}/out/{library}.npz` instead of `{library}/report/fingerprint/{library}.npz`.
    The task also compares them across its libraries by relatedness and var(het), with plots in the report.
    ngs-chew guesses the release from the contig lengths of complete GRCh37, hg19 or GRCh38 references.
    `extra_arguments: ["--genome-release", ...]` sets the release instead, but the reference still needs all contigs of the release.
    So the e2e pipelines with their chr12 reference have no fingerprint task.

22. `target_coverage_report` of `ngs_mapping` becomes an `ngs_data_qc` task with `tool: alfred` and `depends_on.kits`.
    `ngs_mapping` no longer takes `depends_on.kits`.
    The report is `{library}/out/{library}.target_coverage.json.gz` instead of `{library}/report/alfred_qc/{library}.alfred.json.gz`.
    `varfish_export` reads it through `depends_on.qc` and writes its BAM QC only when that field is set.
    Pedigree members from the background libraries have no target coverage there, because `ngs_data_qc` checks no background libraries.
    Main's mapping task did.

### Project section and result names

23. `config.yaml` has a required `project` section, which `snappy init` writes.
    It holds the project's `name` and optionally a `description` and an `external_id`.
    `external_id` is the project's ID in a data management system, such as the SODAR project UUID, which tracks where the raw data lives.

24. The results of `group_by: family`, `donor` and `related` tasks are named after the group instead of the index library or the primary tumor, without a prefix: `work/<family_id>/` or `work/<donor_name>/`.
    This affects joint germline calls, SV calls, their annotation, filtration and checks, and the exports.
    The `write_pedigree.<group>` PED carries the same name in its family column.
    Rerun these tasks, or rename their `work/` directories and files.

### Calling, export and expression

25. `variant_calling` with `tool: bcftools_call` needs `bcftools_call.assembly` (`GRCh37` or `GRCh38`).
    The assembly sets the ploidy of the sex chromosomes for `bcftools call --ploidy`.
    Main guessed it from the path of the reference and called every contig diploid when the path named neither build.

26. `variant_export_external` has no `external_tool` field any more.
    The merged VCF of a pedigree is named after the pedigree, without main's `dragen.` prefix.

27. A `gene_expression_quantification` task with `tool: strandedness` takes a forced strand from `strandedness.strand`, as `ngs_mapping` does with STAR.
    Main took the forced strand from the step's `strand`.
    The step's `strand` applies only to `featurecounts`, `dupradar` and `salmon`.
    It takes 0, 1 or 2 to force unstranded, forward or reverse strand, and -1 to take the decision of `depends_on.strandedness` (`salmon`: to guess).
    A `tool: strandedness` task with a `strand` other than -1 fails when the project loads.

(migrations_example)=

## Example: a neoepitope config from main and its translation

From main's workflow test (`step_config` of the steps involved):

```yaml
static_data_config:
  reference:
    path: /path/to/ref.fa
  cosmic:
    path: /path/to/cosmic.vcf.gz
  dbsnp:
    path: /path/to/dbsnp.vcf.gz
  features:
    path: /path/to/gencode.gtf

step_config:
  ngs_mapping:
    tools:
      dna: ['bwa']
      rna: ['star']
    star:
      path_index: /path/to/star/index
    bwa:
      path_index: /path/to/bwa/index.fa

  gene_expression_quantification:
    tools: [salmon]
    salmon:
      path_index: /path/to/salmon/index

  somatic_variant_calling:
    tools: [mutect2]
    mutect2:
      contamination:
        common_variants: /path/to/common/variants

  hla_typing:
    path_link_in: ''  # OPTIONAL Override data set configuration search paths for FASTQ files
    tools:
      dna: [optitype]   # REQUIRED - available: 'optitype' and 'arcashla'
      rna: [optitype, arcashla]
    optitype:
        max_reads: 5000
        num_mapping_threads: 4

  somatic_variant_annotation:
    path_somatic_variant: ../somatic_variant_calling
    tools: ["vep"]
    is_filtered: false
    vep:
        cache_dir: /path/to/dir/cache

  germline_variant_calling:
    path_ngs_mapping: NGS_MAPPING
    tools: [gatk4_hc]
    gatk4_hc:
      num_threads: 8

  germline_variant_filtration:
    path_variant: ../germline_variant_calling
    has_annotation: false
    filter_list:
      - bcftools:
          exclude: "AD[0:0]+AD[0:1]<50 | AD[0:1]<5 | AD[0:1]/(AD[0:0]+AD[0:1])<0.05"

  combine_variants:
    somatic_variant_type: annotation
    path_somatic_variant: SOMATIC_VARIANT_ANNOTATION
    tool_somatic_variant_annotation: vep
    germline_variant_type: filtration
    path_germline_variant: GERMLINE_VARIANT_FILTRATION
    is_germline_variant_filtered: true
    rename_combined: tumor

  somatic_neoepitope_prediction:
    tools: [pvacseq, pvacfuse, pvacsplice]
    is_filtered: false
    pileup:
      enabled: true
    quantification:
      enabled: true
      duplicate_transcripts_table: /path/to/duplicates.tsv
    phasing:
      enabled: true
      path_combine_variants: COMBINE_VARIANTS
    proteome:
      enabled: true
      path_germline_variants: GERMLINE_VARIANT_FILTRATION
      germline_variant_step: germline_variant_filtration
      external_proteome: /path/to/gencode.fa
    tools_hla_typing:
      dna:
        class_i: [optitype]
      rna:
        class_i: [optitype]
        class_ii: [arcashla]
    pvacseq:
        class_ii_epitope_length: [10, 11]
        extra_args: ["--percentile-threshold-strategy exploratory"]
        net_chop:
          enabled: true
          path_netchop: /path/to/netchop.bin
        netmhc_stab:
          enabled: true
    pvacfuse:
        algorithms: [all_class_i]
    pvacsplice:
        algorithms: [NetMHCIIpan, MHCnuggetsII]
        use_all_transcripts: true
        genes_of_interest_file: /path/to/genes.txt
        net_chop:
          enabled: true
          path_netchop: /path/to/netchop.bin
        netmhc_stab:
          enabled: true
        
data_sets:
  first_batch:
    file: sheet.tsv
    search_patterns:
    - {'left': '*/*/*_R1.fastq.gz', 'right': '*/*/*_R2.fastq.gz'}
    search_paths: ['/path']
    type: matched_cancer
    naming_scheme: only_secondary_id
```

The same config in snake9 loads with `load_project()`.
On a test sheet with a tumor RNA library, its dry-run builds 84 jobs with either phasing tool:

```yaml
project:
  name: neoepitope

sample_sheets: [sheet.tsv]

tasks:
  # The FASTQ files of the libraries, listed by snappy scan-files
  - name: reads
    step: external_data
    config:
      produces: {type: reads}
      library_files: [reads.tsv]

  # Reference data (static_data_config and the path_* fields on main)
  - name: genome
    step: external_data
    config:
      produces: {type: sequences, tags: [reference, dna]}
      files: {fasta: /path/to/ref.fa}
  - name: genes
    step: external_data
    config:
      produces: {type: features}
      files: {gtf: /path/to/gencode.gtf}
  - name: dbsnp
    step: external_data
    config:
      produces: {type: variants, tags: [dbsnp]}
      files: {vcf: /path/to/dbsnp.vcf.gz}
  - name: exons
    step: external_data
    config:
      produces: {type: regions, tags: [exons]}
      files: {bed: /path/to/exons.bed}
  - name: bwa_index
    step: external_data
    config:
      produces: {type: index, tags: [bwa, dna]}
      files: {index: /path/to/bwa/index.fa}
  - name: star_index
    step: external_data
    config:
      produces: {type: index, tags: [star, rna]}
      files: {index: /path/to/star/index}
  - name: salmon_index
    step: external_data
    config:
      produces: {type: index, tags: [salmon, rna]}
      files: {index: /path/to/salmon/index, transcript_to_gene: /path/to/tx2gene.tsv}
  - name: common_variants
    step: external_data
    config:
      produces: {type: variants, tags: [common]}
      files: {vcf: /path/to/common/variants.vcf.gz}
  - name: vep_cache
    step: external_data
    config:
      produces: {type: database, tags: [vep]}
      files: {cache: /path/to/dir/cache}
  - name: proteome
    step: external_data
    config:
      produces: {type: sequences, tags: [proteome]}
      files: {fasta: /path/to/gencode.fa}

  # Mapping and expression
  - name: mapping_dna
    step: ngs_mapping
    config:
      depends_on: {reads: reads, index: bwa_index, reference: genome}
      tool: bwa
      bwa: {}
  - name: mapping_rna
    step: ngs_mapping
    config:
      depends_on: {reads: reads, index: star_index, reference: genome, features: genes, exons: exons}
      tool: star
      star: {}
  - name: expression
    step: gene_expression_quantification
    config:
      depends_on: {reads: reads, features: genes, index: salmon_index}
      tool: salmon
      salmon: {}
  - name: fusions
    step: somatic_gene_fusion_calling
    config:
      depends_on: {reads: reads, index: star_index, reference: genome, features: genes}
      tool: arriba
      arriba: {}

  # HLA typing: one task per tool and molecule (tools_hla_typing on main)
  - name: hla_optitype_dna
    step: hla_typing
    config:
      depends_on: {reads: reads}
      tool: optitype
      library_selection: "extraction_type == 'dna'"  # tumor and normal
      optitype: {max_reads: 5000, num_mapping_threads: 4}
  - name: hla_optitype_rna
    step: hla_typing
    config:
      depends_on: {reads: reads}
      tool: optitype
      library_selection: "extraction_type == 'rna'"
      optitype: {max_reads: 5000, num_mapping_threads: 4}
  - name: hla_arcashla
    step: hla_typing
    config:
      depends_on: {alignments: mapping_rna}
      tool: arcashla
      arcashla: {}

  # Somatic calls, annotated by VEP with the pVACtools plugins
  - name: calling
    step: variant_calling
    config:
      depends_on: {alignments: mapping_dna, reference: genome, common_variants: common_variants}
      tool: mutect2
      mutect2:
        contamination: {enabled: true}
  - name: annotation
    step: variant_annotation
    config:
      depends_on: {variants: calling, reference: genome, database: vep_cache}
      tool: vep
      vep: &vep_pvactools
        plugins:
          - {name: Frameshift, path: /path/to/vep_plugins/Frameshift.pm}
          - {name: Wildtype, path: /path/to/vep_plugins/Wildtype.pm}
  - name: annotation_all_transcripts  # use_all_transcripts on main
    step: variant_annotation
    config:
      depends_on: {variants: calling, reference: genome, database: vep_cache}
      tool: vep
      vep: &vep_all_transcripts
        <<: *vep_pvactools
        pick: false

  # Germline calls of the normals, combined with the somatic calls and phased with the tumor DNA
  - name: germline_calling
    step: variant_calling
    config:
      depends_on: {alignments: mapping_dna, reference: genome, dbsnp: dbsnp}
      tool: gatk4_hc_joint
      library_selection: &normals "role == 'normal' and extraction_type == 'dna'"
      group_by: library
      gatk4_hc_joint: {num_threads: 8}
  - name: germline_annotation
    step: variant_annotation
    config:
      depends_on: {variants: germline_calling, reference: genome, database: vep_cache}
      tool: vep
      library_selection: *normals
      vep: *vep_all_transcripts
  - name: germline_filtration
    step: variant_filtration
    config:
      depends_on: {variants: germline_annotation}
      tool: bcftools
      library_selection: *normals
      bcftools:
        exclude: "AD[0:0]+AD[0:1]<50 | AD[0:1]<5 | AD[0:1]/(AD[0:0]+AD[0:1])<0.05"
  - name: combined
    step: combine_variants
    config:
      # All transcripts: pVACseq matches the nearby variants by the somatic variant's transcript
      depends_on: {somatic_variants: annotation_all_transcripts, germline_variants: germline_filtration, reference: genome}
      rename_combined: tumor
  - name: phasing
    step: variant_phasing
    config:
      depends_on: {combined_variants: combined, alignments: mapping_dna, reference: genome}
      tool: gatk_read_backed_phasing_tumor  # as on main; or whatshap_tumor

  # Neoepitopes: pvacsplice reads all transcripts, so it is a task of its own
  - name: neoepitopes
    step: somatic_neoepitope_prediction
    config:
      depends_on:
        variants: annotation
        hla_types: &hla_types [hla_optitype_dna, hla_optitype_rna, hla_arcashla]
        reference: genome
        alignments: mapping_rna
        phased_variants: phasing
        expression: expression
        fusions: fusions
        proteome: proteome
      tool: pvactools
      pvactools:
        pvacseq:
          class_ii_epitope_length: [10, 11]
          extra_arguments: ["--percentile-threshold-strategy exploratory"]
          pileup: {enabled: true}
          expression: {enabled: true, duplicate_transcripts_table: /path/to/duplicates.tsv}
          phasing: {enabled: true}
          net_chop: {enabled: true, path: /path/to/netchop.bin}
          netmhc_stab: {enabled: true}
        pvacfuse:
          enabled: true
          algorithms: [all_class_i]
  - name: neoepitopes_splice
    step: somatic_neoepitope_prediction
    config:
      depends_on:
        variants: annotation_all_transcripts
        hla_types: *hla_types
        reference: genome
        features: genes
        alignments: mapping_rna
        proteome: proteome
      tool: pvactools
      pvactools:
        pvacseq: {enabled: false}
        pvacsplice:
          enabled: true
          algorithms: [NetMHCIIpan, MHCnuggetsII]
          genes_of_interest: /path/to/genes.txt
          net_chop: {enabled: true, path: /path/to/netchop.bin}
```

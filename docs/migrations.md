# Migrating from snappy 0.6

This page collects what changed between snappy 0.6 (the `main` branch) and the task-based snappy on Snakemake 9: the configuration model in detail, the new command and defaults, what was removed or changes results on purpose, and step-by-step migration notes for existing projects.
The pipeline configs in `.tests/test-workflow/pipelines/*/config.yaml` are complete examples of the new layout.

## Task-based configuration

Instead of one block per workflow step, `config.yaml` lists **named tasks**.
Each task names a step type (`ngs_mapping`, `variant_calling`, ...) and carries its own config: the tool, the library selection, the dependencies and the tool settings.
Tasks form a DAG through their `depends_on` fields, and `external_data` tasks provide existing files such as the genome, enrichment kit regions, indexes and databases.

By default `snappy run` builds every task and links its result files into `output/`.
`snappy run --task <name>` builds one task, and `--frozen` with `--task` uses the existing outputs of the upstream tasks.

## Single tool per task

Workflows used to register all available tool implementations at once.
Now each workflow has a single `tool` field, an enum that selects exactly one implementation.
Tool-specific configuration is nested under the tool name directly:

```yaml
config:
  tool: gatk4_hc_gvcf
  gatk4_hc_gvcf:
    window_length: 10000000
```

This replaces the old multi-tool list pattern.
It simplifies the Snakefiles, makes config validation more precise, and avoids ambiguity about which tool's resources or output paths apply.
If a workflow needs several tools that cannot be split into separate tasks (rare, but possible), the workflow code handles that; the single-tool constraint applies at the config and task level.

## Library selection

Every task can filter which libraries it processes with a pandas `DataFrame.query()` expression:

```yaml
library_selection: "role == 'tumor' and extraction_type == 'dna'"
```

The expression runs against a DataFrame with one row per NGS library.
The columns that are **always available**:

| Column            | Description                                                                    |
| ----------------- | ------------------------------------------------------------------------------ |
| `library_name`    | NGS library identifier                                                         |
| `extraction_type` | "dna", "rna", ... (always lowercased)                                          |
| `kind`            | "cancer" or "germline"                                                         |
| `role`            | "tumor" / "normal" / "index" / "father" / "mother" / "affected" / "unaffected" |
| `is_primary`      | True for the primary tumor (cancer) or the pedigree index (germline)           |
| `sex`             | "male", "female", or "unknown"                                                 |
| `donor_name`      | Patient or family identifier                                                   |
| `sample_name`     | Bio-sample name                                                                |
| `cohort_name`     | Donor name for cancer, pedigree name for germline                              |
| `tissue_type`     | "tumor", "normal", or "unknown"                                                |
| `disease_state`   | Integer disease state (0 = unaffected, 1 = affected)                           |
| `father_name`     | Parent library name (or "0")                                                   |
| `mother_name`     | Parent library name (or "0")                                                   |
| `data_set`        | Name of the library's data set in the config                                   |
| `is_background`   | True for libraries of background data sets                                     |

Further sample sheet columns are added in snake_case, for example `libraryKit` as `library_kit` and `libraryType` as `library_type`.

Columns that are **optionally available** (via the `relationships` config, which merges with per-step defaults):

- Any column added by a `RelationshipDefinition` in the step's `default_relationships` or the task's config-level `relationships`.

A task's libraries come from three layers:

- What its step can process at all, which always applies: for example RNA libraries for expression quantification and fusion calling, tumor DNA libraries for mutect2, normal DNA libraries for a panel of normals, and only libraries of cancer sample sheets for the somatic steps.
- Without `library_selection`, the step's default: `role == 'tumor' and extraction_type == 'dna'` on cancer sample sheets and `extraction_type == 'dna'` on germline sample sheets, unless the step declares its own (mapping takes all libraries, the RNA steps all RNA libraries, a panel of normals all normal DNA libraries).
- `library_selection` replaces that default, for example to select a batch (`data_set == 'batch2'`) or to drop samples, but never goes beyond what the step can process.

Several data sets may have libraries of the same name; a task that would take two of them fails and suggests a `library_selection` such as `data_set == 'batch1'`.

Libraries of data sets with `is_background: true` (a reference cohort, for example) are no task's results.
Mapping and panels of normals process them, and the SV and gCNV callers add them as cohort members or model samples; `background_selection` narrows them per task, such as `data_set == 'cohort_1000g'`.

(Note: the DataFrame is built from biomedsheets shortcut objects.
We plan to move away from biomedsheets as a dependency in a future release, but the DataFrame interface and the query-based selection mechanism will stay.)

## Group-by for aggregation level

Some steps (variant calling, annotation, filtration) can produce output at different levels of aggregation:

- `group_by: cohort`: one output per pedigree / family; for a germline trio, one joint VCF for the whole family, named after the index library.
  For a cancer project, one joint VCF per donor (tumor and normal together).
- `group_by: library`: one output per library.
- Unset: the step's default, else the setting that upstream tasks declare, else `library`.
  Germline variant calling defaults to `cohort` and calls each pedigree jointly.

The setting propagates **downstream**: a filtration step after cohort-wise variant calling writes per cohort too, unless it sets `group_by` itself.
This keeps the aggregation level consistent through a chain of dependent tasks.
Steps that read one VCF per pedigree (VarFish export, peddy, IGV sessions, pedigree phasing) check when the project loads that their upstream task writes per cohort.

## Dependencies and task-local files

A task gets its input files in one of two ways.
A **dependency** comes from another task, named in `depends_on`.
It has a data type and tags, and the consumer checks both, and the keys it reads, when the project loads.
Reference data are dependencies too: the genome, the gene annotation, dbSNP, enrichment kits, exons and other regions, and the mapper indexes.
`external_data` tasks provide existing files, and `reference_index` tasks build indexes.
So each reference file is named once, and a switch to another genome build only changes these tasks.

A file that only one tool reads can still be a dependency.
A bwa index has a format that only bwa reads, but it must match the genome that the other tasks use.
So it is `depends_on.index`, from a task with `produces: {type: index, tags: [bwa, dna]}`.

A **task-local file** belongs to the setup of one tool and does not depend on the reference, such as a jar or a Python file with vembrane helpers.
It is a path in the tool section, which the step passes to the rule as an input.

BED files of regions have their own data type `regions` with the key `bed`.
The regions filter accepts any of them in `depends_on.regions`, and steps that need a certain kind of regions require a tag (`exons`, `kits`, `mappable`, `microsatellites`, ...).
Population variants are `variants` with a tag, like dbSNP: `depends_on.germline_resource` (tag `population`) and `depends_on.common_variants` (tag `common`).
Every tool that reads a GTF reads `depends_on.features`; RNA-SeQC needs a GTF prepared for it, so it requires the tag `rnaseqc` as well.
A tool's own database comes from `depends_on.database`, a task of type `database` tagged with the tool, such as the VEP cache or the mehari databases.
The field's requirement names the task's tool, which the loader fills in, so a VEP cache wired to a mehari task fails when the project loads.
One mehari database task can serve both variant annotation and the VarFish export.
Every job that reads the genome also has its `.fai` and `.dict` as inputs, so a missing sidecar fails when Snakemake builds the DAG.
Likewise, a job that reads dbSNP or population variants has the VCF's index as an input: `.tbi` for a `.vcf.gz`, `.idx` for a plain `.vcf`.

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

All reference-specific files of the maintained tools follow this rule now.
The external export steps (`variant_export_external`, `wgs_cnv_export_external`, `wgs_sv_export_external`) are left for later: they still run varfish-annotator with Jannovar `.ser` files from path fields.

## Relative path resolution

Every path in `config.yaml` is relative to the project directory, the directory of the config file: sample sheets, search paths, `external_data` files and files in tool sections.
`snappy run` passes the config as `--configfile`, so this also holds when Snakemake runs in another directory.
The top-level `pathvars` section names roots for paths, such as `<refs>/genome.fa`, with the defaults and semantics of Snakemake's `pathvars`.
snappy does not check whether the files exist when the project loads, so dry-runs work without them; `snappy check-files` lists the missing ones.

## Work and output directories

Each task has its own directory `tasks/<task>/`.
`work/` holds every file the jobs write and is the contract between tasks: a task reads the `work/` files of its upstream tasks, never their `output/`.
`output/` holds only the result files: one local `<task>_results` job per task with result files links them from `work/` and writes `output/.manifest.tsv`.
Only result files get an `.md5` file, logs hold only the log, and `snappy run` archives the logs of each run.
Scatter fragments are temporary, and `delete_intermediates: [intermediates]` makes per-sample intermediates temporary, too; `[bam]` does so for the BAM of a mapping task that a later task processes further.

## Developer contracts: produces and DataSignature

For developers writing new workflow steps, `load_project()` validates the dependency chains before Snakemake builds any rule.

- **`produces`** (or the classmethod `task_produces(config, upstream)`) declares the `DataSignature` values of a task's outputs (e.g. `ALIGNMENTS` with tag `dna`).
- **`depends_on` fields** carry the `DataSignature` they require and an `ExpectedPathSchema` of the files they read.
- **`DataSignature`**: a frozen descriptor pairing a `DataType` enum with a frozenset of tags.
  Tag constraints support string literals, "OR" tuples, and "NOT" prefixes (e.g. `-"somatic"`).

`load_project()` checks that each upstream task produces what the `depends_on` field requires, and that an `external_data` task has the keys the consumer reads.
This catches mismatches like wiring an RNA alignment step into a variant caller that expects DNA, or a germline caller into a somatic pipeline.
It also checks that an index from a `reference_index` task is built from the genome the reader names, and that `FILTER` IDs cannot mean two things.
Consumers read upstream files only through `get_upstream_paths()`, which replaces the old implicit "the paths will probably work out" approach.

## The `snappy` command

The old entry-point scripts are replaced by a single `snappy` command with subcommands:

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

### `snappy run` options

- `--task <name>`: build only the named task and the upstream files it needs, and link its result files (default: every task).
- `--frozen`: with `--task`, use the existing outputs of the upstream tasks and never rebuild them.
- `--directory`: the project directory (default: the current directory).
- `--verbose`: more log output.
- Everything after `--` is passed verbatim to Snakemake.

## Changed defaults

### Default workflow profile

A default Snakemake profile is always applied.
It sets:

- `software-deployment-method: conda`: conda environments are enabled by default, no need for `--use-conda`.
- `logger: snkmt`: structured logging to an SQLite database at `.snakemake/log/snkmt.sqlite`.
  Watch it live with `snappy watch`.
- `persistence-backend: db`: a single SQLite file for Snakemake's metadata instead of many small files in `.snakemake/`.
  This matters when running on shared network filesystems (CephFS, NFS) that handle a single large file much better than thousands of tiny ones.

### Cluster execution

The executor comes from a Snakemake profile of the system or the user (`--profile <name>` after `--`, or `SNAKEMAKE_PROFILE`), not from snappy: the old `--profile-snappy-pipeline` flag, the `pipeline_job.sh` job script and the `snappy refresh` and `snappy status` commands are gone.
The Snakemake process runs in a terminal multiplexer on the login node or in an allocation (`salloc`), as the documentation of the SLURM executor recommends; {ref}`usage` has an example profile for SLURM.

### Output path patterns

Result file paths are now much shorter.
The old patterns embedded the tool name and step substage (e.g. `bwa.filtered.vcf`, `bwa.raw.vcf`) into every path, even when there was only one tool and one substage involved.
Since those values are already defined by the task config, repeating them in paths was redundant.
Paths now contain only the truly variable parts, typically just the library name (or cohort name).
In multi-tool workflows or substeps that require disambiguation, additional wildcards are still used, but the default case is a single wildcard.

Profiles are now passed via `--workflow-profile`, which is the recommended approach in Snakemake 9.
Multiple `--workflow-profile` arguments are layered (later profiles override earlier ones), so the default profile is applied first and the SLURM profile overrides only what it needs to.

## Snakemake 7 to Snakemake 9

The pipeline now requires Snakemake >= 9 and no longer supports Snakemake 7.
The most visible change is that **subworkflows** have been replaced by Snakemake 9's **module** system.
Each step's Snakefile is imported as a module with a `replace_prefix` mapping, so output and work directories are namespaced per task (e.g. `tasks/sample_bwa/output/` instead of a flat `output/`).
This is what makes it possible to run the same step type multiple times in one project without file collisions.

More wrappers were migrated to the typed wrapper classes (`RWrapper`, `PythonWrapper`, `ShellWrapper`).
These improve type checking and make wrapper signatures more explicit compared to the plain `Wrapper` base class.

The Snakemake version bump also required updates to several dependencies (pandas, pydantic, etc.) to stay compatible with Snakemake's own dependency chain.

snappy uses the development version of Snakemake, which deploys software through plugins.
The conda plugin builds the wrapper environments from their pin files, below `.snakemake/software/`; `snappy run -- --deploy-software-envs` creates them without running jobs (formerly `--conda-create-envs-only`).
R packages that are not on conda-forge or bioconda in the needed versions (for sequenza and scarHRD) are installed by the post-deploy script of the wrapper's environment, at fixed versions or commits, when the environment is created.
PureCN and pVACtools run from their container images (`software: container(...)`), which Snakemake pulls and runs with apptainer.

## Project management with pixi

We adopted [pixi](https://pixi.sh) as the project manager (install pixi: see <https://pixi.sh/latest/#installation>).
Snappy has both conda dependencies (for system tools like BWA, STAR, samtools) and PyPI dependencies (for Python packages not available on conda), and pixi handles both ecosystems natively.

The pixi configuration defines several development tasks:

- `pixi run -e dev test`: run the test suite
- `pixi run -e dev lint`: run all linters (ruff check, ruff format --check, snakefmt --check).
  Note: `mypy` is not currently in use.
- `pixi run -e dev fmt`: auto-format with ruff
- `pixi run -e dev snakefmt`: auto-format Snakemake files
- `pixi run -e dev srcfmt`: run both formatters (depends-on: fmt, snakefmt)

Pixi also replaces the old `conda env create` / `pip install` workflow for developers, a single `pixi install` sets up everything.

## Build system and packaging

The build backend was migrated from `setuptools` to `hatchling`.
The old `setup.py`, `setup.cfg`, `environment.yml`, and `Makefile` have been removed; pixi handles environment management and the build configuration lives entirely in `pyproject.toml`.

## CI and test suite modernization

The CI pipeline was updated to use `prefix-dev/setup-pixi` instead of the legacy `conda-cache` GitHub action.
End-to-end integration tests (`.github/workflows/ci-e2e.yml`) now run test workflow configurations directly with `pixi run snappy run` against dedicated test datasets checked in via Git LFS (`.tests/test-workflow/pipelines/`).
They dry-run and then run the cancer, germline and reference pipelines on small test data.
Wrapper environments are pinned (`environment.linux-64.pin.txt`); in a pull request, a changed environment file must come with its changed pin file (`pixi run -e dev update-pins <file>`).

## What was removed

- **somatic_wgs_sv_calling workflow**: this workflow was a specialised SV calling pipeline for whole-genome data with its own complex config.
  It was not widely used and has been removed entirely.
  SV calling for targeted data via `sv_calling_targeted` is unaffected.
- **Old `tpls/profile/` directory**: superseded by the workflow profile `snappy_pipeline/apps/profile/` and a Snakemake profile of the system or the user for the executor.
- **`--profile-snappy-pipeline` flag, `pipeline_job.sh`, `snappy refresh` and `snappy status`**: cluster execution is configured in a Snakemake profile.
- **Multi-tool list configs**: the `tools:` list in step configs is no longer supported.
  Use `tool:` with a single value.
- **Legacy `snappy-snake` entry point**: removed; use the `snappy` command instead.
- **R package install jobs of sequenza and scarHRD** and their `work/R_packages/` directories: the environments install these packages when they are created.
- **Container builds of PureCN and pVACtools**: the rules `panel_of_normals_purecn_install` (the `singularity` wrapper, `work/containers/`) and `somatic_neoepitope_prediction_pvactools_install` (`work/container/`), and the key `container` of a PureCN panel's outputs; Snakemake pulls the images.
- **Jannovar statistics report** (`jannovar_stats` of `variant_calling`).
- **`mbcs` tool of `ngs_mapping`** and its nested Snakemake run (`nested_snakemake.py`), with `depends_on.common_variants` of `ngs_mapping`: a chain of tasks replaces it (migration note 20).
- **`fusioncatcher` tool of `somatic_gene_fusion_calling`**: it had no wrapper (main neither), so it could not run.
- **`defuse`, `hera`, `jaffa` and `pizzly` tools of `somatic_gene_fusion_calling`**: untested and probably broken; `arriba` (recommended) and `star_fusion` remain.
- **main's germline filters** (`variant_denovo_filtration`, the germline filters of main's `variant_filtration` and the `protected` rescue filter): they no longer ran on main (a deleted shell library, Jannovar inputs); `vembrane` expressions in a `variant_filtration` task cover their per-record criteria.
- **Config fields that nothing read or that duplicated the upstream config**: `cbioportal_export` `vcf2maf.annotation_tool` (the export requires VEP annotations), `copy_number_alteration.copy_number_tool` and `patient_info`; `somatic_cnv_checking` `cnv_assay_type`; `hla_typing` `arcashla.mapper`; `somatic_targeted_seq_cnv_calling` `somatic_variant_caller` of PureCN.
  A config that sets one fails validation.

## What did not change, and deliberate result changes

The rules run the same tools (BWA, STAR, GATK, etc.).
A real-run comparison of the e2e pipelines against main 0.6.3 gives the same results, with two exceptions: main's bcftools 1.15.1 drops every record of the dkfz-filtered calls (it cannot parse their GT), while snake9's bcftools 1.23.1 keeps and tags them, and `FILTER` IDs are task names (migration note 10).
The per-step Snakefiles (`workflows/*/Snakefile`) still use the same patterns: `wf.config`, `wf.get_input_files()`, `wf.get_params()` and wrappers.
Outside the e2e pipelines, results change on purpose where main was wrong or ignored a setting: MELT's GroupAnalysis reads the discovery output of all pedigree members (main read the first member's only), the `ignore_contigs` defaults apply to mutect2 and sequenza (migration note 14), arcasHLA gets `zygosity_threshold` (note 17), VEP annotates multiallelic records split into one record per allele (`bcftools norm`), the ebfilter panel draws from all normal DNA libraries instead of the normals matched to a tumor, Control-FREEC takes each donor's sex from the sample sheet (main used `XX` for every sample), a PureCN panel reads the GenomicsDB of its `mutect2` panel (main fails there, because the `genomicsdb` R package of the PureCN image needs integers where GATK writes strings), the BWA mappers read a fixed 100 Mbp per batch (`bases_per_batch`, `-K`), so their alignments no longer depend on `num_threads_align`, their read group IDs are the names of the discovered read groups (main: `<library>.0`, `<library>.1`, ... for `bwa` and one read group for all reads for `bwa_mem2`), and `bwa_mem2` maps single-end reads (main read the right file, which they do not have) and fails when it cannot decompress the reads (main wrote an empty BAM).

## Migration notes

To migrate an existing project:

1. Wrap each step's configuration in a `tasks:` list entry with `name` and `step` fields.
2. Replace `tools: [bwa]` with `tool: bwa` and move tool-specific config out of the list to a top-level key matching the tool name.
3. Remove `--profile-snappy-pipeline` references from scripts; the workflow profile of snappy is applied automatically.
   Put the cluster settings into a Snakemake profile of the system or the user and pass `--profile <name>` after `--`, or set `SNAKEMAKE_PROFILE`.
4. Use `snappy run` instead of the old entry points.
5. The old `config.yaml` structure without a `tasks:` key is rejected; run `snappy init` in a fresh directory to see the new layout.
6. Move files that depend on the reference out of tool sections into `external_data` tasks named in `depends_on` (see Dependencies and task-local files); for example, a VEP task now needs a cache task in `depends_on.database`.
7. Sample sheets are no longer looked up in the parent directory of the project, and `data_sets` search paths are relative to the project directory instead of the sample sheet's directory.
   Give such paths relative to `config.yaml`, absolute, or with a `pathvars` root.
8. `data_sets` `search_patterns` are regular expressions for the path below a library's folder, with a `(?P<readgroup>...)` group that pairs the `left` and `right` files, instead of glob patterns.
   A data set has no default `search_paths` (`../raw`) any more; `search_patterns` default to `(?P<readgroup>.+)\.R1\.fastq\.gz` and `(?P<readgroup>.+)\.R2\.fastq\.gz`.
   A library's reads are searched only below the search paths of its own data set, not those of all data sets.
9. VEP annotates all transcripts in a task of its own with `vep.pick: false`, instead of the `.full.vcf.gz` side output; consumers that need all transcripts (`use_all_transcripts` of neoepitope prediction, `use_all_transcripts_for_*` of combine_variants) read such a task.
   VEP plugins (`vep.plugins`) are local `.pm` files; main's download from a URL is gone.
   Annotations are tagged with the annotator (`vep` or `mehari`), so a consumer can require one.
10. The `bcftools`, `regions` and `ebfilter` filters of `variant_filtration` use the task name as FILTER ID instead of main's `<filter>_<position in filter_list>` (`bcftools_2`); expressions that test FILTER values need the new names.
    No task may be named `PASS` or like a FILTER ID that another task writes (the filters of mutect2, `LowQual` of the GATK callers, the expression names of `vembrane` in tag mode, `bPcr` and `bSeq` of `dkfz`), and no task may write an ID that a task upstream of it writes.
11. Every tool's config section takes `extra_arguments`, a list of additional arguments for the tool's main command line; a single-tool step without a section (such as `variant_checking` or `repeat_expansion`) takes it on its step model.
    The former `extra_args` fields (`bwa`, `bwa_mem2`, AGeNT, pVACtools, `reference_index` and `vembrane`, whose value was a string) are renamed.
    A tool with several commands passes them to one, such as HaplotypeCaller for `gatk4_hc_gvcf`, `segment` for cnvkit or `GermlineCNVCaller` for gCNV.
    The tumor phasing tools (`gatk_read_backed_phasing_tumor`, `whatshap_tumor`) and `somatic_variant_signatures` (`tabulate_vcf`) have new sections for it.
12. The BED file of mappable regions for gCNV needs its index next to it (`.idx` from `gatk IndexFeatureFile`, or `.tbi` for a bgzipped BED); `snappy check-files` lists it.
13. `somatic_neoepitope_prediction` is one task with `tool: pvactools`; `pvactools.pvacseq` (on by default), `pvacfuse` and `pvacsplice` are toggles instead of the `tools` list.
    Its inputs come from `depends_on` instead of `path_*` and `tool_*` fields: `hla_types` lists the HLA typing tasks instead of `tools_hla_typing`, `alignments` (RNA) replaces the mapper of `pileup`, `pileup`, `quantification` (now `expression`) and `phasing` are pvacseq options, and `path_netchop` is `net_chop.path`.
    Phasing runs in a `variant_phasing` task with `tool: gatk_read_backed_phasing_tumor` (GATK 3 ReadBackedPhasing of the combined variants with the tumor DNA, as on main) or `whatshap_tumor` (WhatsHap, writing ReadBackedPhasing's `HP` tag, which pVACseq reads), and the neoepitope task reads it through `depends_on.phased_variants`.
    The phased variants should be annotated for all transcripts, as pVACseq matches nearby variants by the transcript of the somatic variant.
    The variants must come from VEP with the pVACtools plugins `Frameshift` and `Wildtype`, which the task checks when Snakemake builds the workflow; for phasing, the germline calls must be VEP-annotated, too.
    Predictions that read different annotations (pvacsplice with all transcripts, pvacseq without) are separate tasks.
    The proteome is an external FASTA (`depends_on.proteome`); building it from germline variants is not ported.
   The example at the end of this page translates a neoepitope config from main.
14. `ignore_chroms` is `ignore_contigs` in every step (`variant_calling`, `variant_phasing`, `sv_calling_wgs`, `panel_of_normals` and the `sequenza` section).
    The entries are regular expressions matched with `re.search`, where main's somatic steps used globs: `HLA-*` becomes `^HLA-`, `*_decoy` becomes `_decoy$`.
    In `variant_calling`, the list is a field of each caller's section instead of the step, so a YAML anchor shares one list between tasks; mutect2 skips decoys, HLA, EBV and the hs37d5 extras by default, like the germline callers (on main it called on every contig).
15. `bbduk.adapter_sequences` takes reference files only, with paths relative to the config file; BBDuk's keywords (`adapters`, `artifacts`, `phix`, `lambda`, `pjet`, `mtst`, `kapa`) go in `bbduk.builtin_references`.
16. The `tag_mode` of main's vembrane filter is gone: an expression names the variants to tag, so `tag_mode: include` becomes the negated expression, `not (...)`.
17. `arcashla.zygocity_threshold` is `arcashla.zygosity_threshold`, and it now reaches arcasHLA (main never passed it, so arcasHLA ran with its own default 0.1 instead of the configured 0.15).
18. `somatic_wgs_cnv_calling` with `tool: cnvkit` needs `depends_on.panel_of_normals`: a whole-genome `panel_of_normals` task with `tool: cnvkit` (a flat reference when there are no normals).
19. `purecn.container` and `pvactools.container` are the images that Snakemake pulls, pinned by digest by default (`markusriester/purecn` of 2024-05-28, `griffithlab/pvactools:7.1.3`); the path of a local `.sif` file works too, and `pvactools.container: null` is no longer valid.
    A PureCN calling task runs the image of its panel of normals task.
    Jobs in a container see the project, `$HOME`, `/tmp` and the snappy installation; other directories they read go into `APPTAINER_BIND`, see {ref}`usage`.
20. `tool: mbcs` of `ngs_mapping` becomes a chain of tasks, see the recipe in {ref}`somatic_ngs`: `adapter_trimming` with `tool: agent` (from `agent.prepare`, `use_barcodes`), a mapping task with the tool of `mbcs.mapping_tool`, `alignment_processing` with `tool: agent` (from `agent.mark_duplicates`) and `alignment_processing` with `tool: gatk4_bqsr` (`recalibrate`, `depends_on.common_variants`).
    `lib_prep_type` and `consensus_mode` are required now; main's mbcs passed `-None` and `-c=None` without them, and it could not run with `bwa` at all (literal tabs in its read group line, which bwa rejects) nor with `recalibrate: false`.
    The read groups are those of the data set's search pattern, as in mbcs (`<library>_L001`); duplicates are masked across the lanes of a library, while mbcs masked them per lane, so more reads carry the duplicate flag.
    With main's mbcs output without BQSR as input, `gatk4_bqsr` writes main's recalibrated BAM.

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

The same in snake9 (loads with `load_project()`; the dry-run builds 84 jobs on a test sheet with a tumor RNA library, with either phasing tool):

```yaml
data_sets:
  cancer:
    file: sheet.tsv
    type: matched_cancer
    search_paths: [raw]

tasks:
  # Reference data (static_data_config and the path_* fields on main)
  - name: genome
    step: external_data
    config:
      produces: {type: raw, tags: [reference, dna]}
      files: {fasta: /path/to/ref.fa}
  - name: genes
    step: external_data
    config:
      produces: {type: raw, tags: [features]}
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
      produces: {type: raw, tags: [proteome]}
      files: {fasta: /path/to/gencode.fa}

  # Mapping and expression
  - name: mapping_dna
    step: ngs_mapping
    config:
      depends_on: {reads: data_sets, index: bwa_index, reference: genome}
      tool: bwa
      bwa: {}
  - name: mapping_rna
    step: ngs_mapping
    config:
      depends_on: {reads: data_sets, index: star_index, reference: genome, features: genes, exons: exons}
      tool: star
      star: {}
  - name: expression
    step: gene_expression_quantification
    config:
      depends_on: {reads: data_sets, features: genes, index: salmon_index}
      tool: salmon
      salmon: {}
  - name: fusions
    step: somatic_gene_fusion_calling
    config:
      depends_on: {reads: data_sets, index: star_index, reference: genome, features: genes}
      tool: arriba
      arriba: {}

  # HLA typing: one task per tool and molecule (tools_hla_typing on main)
  - name: hla_optitype_dna
    step: hla_typing
    config:
      depends_on: {reads: data_sets}
      tool: optitype
      library_selection: "extraction_type == 'dna'"  # tumor and normal
      optitype: {max_reads: 5000, num_mapping_threads: 4}
  - name: hla_optitype_rna
    step: hla_typing
    config:
      depends_on: {reads: data_sets}
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

(usage)=

# Usage

The `snappy` command is the interface to the pipeline.
`snappy --help` lists the subcommands.

(usage_commands)=

## Subcommands

The table is generated from the help of each command, which `snappy <command> --help` prints in full.

```{snappy-commands}
```

`snappy init` and `snappy task add` take `[task_name=]step_type` arguments, for example `mapping=ngs_mapping`.
Besides `config.yaml`, `snappy init`:

- writes an empty sample sheet (`--sheet-type generic|germline|cancer`) that `config.yaml` names,
- records the project name, external id and description,
- makes the project a git repository that tracks only `.gitignore`, `config.yaml` and `README.md` (`--track-sheet` adds the sheet; `--remote <url>` adds `origin` without pushing; `--no-git` skips it),
- asks on a terminal for what no option gives (`--no-interactive` takes the defaults).

`--config-mode` sets the style of the generated task configuration: `commented` (default; every field with its description, the optional ones commented out), `full` (every field active) or `minimal` (the required fields only, without descriptions).
`config.yaml` lists `pathvars`, `data_sets` and `resources` before `tasks`, so that added tasks end the file; with `minimal`, the header of `config.yaml` points to the docs instead of explaining the sections.

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
: Also log debug messages, such as the Snakemake command line.

Everything after `--` goes to Snakemake:

```shell
snappy run -- --dry-run --cores 1               # list the jobs, run nothing
snappy run --task mapping -- --cores 8          # one task and its missing upstream files
snappy run --task calling --frozen -- --cores 8 # reuse the existing mapping, never rerun it
snappy run -- --profile slurm                   # all tasks through the executor of a profile
```

## Cluster execution

`snappy run` always uses the workflow profile of snappy.
It deploys the software with conda and, for the tools that run from a container image, with apptainer, and sets the `snkmt` logger and the database persistence backend.

The executor and its settings belong to the system or the user, in a Snakemake profile:

- below `/etc/xdg/snakemake/` for the system, or `~/.config/snakemake/` for the user,
- selected with `--profile <name>` after `--`, or with the environment variable `SNAKEMAKE_PROFILE`.

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

### Job groups

On a cluster, every job is a cluster job, unless its rule is local.
The task field `job_group` puts the jobs of a task into a Snakemake job group.
Dependent jobs of the tasks with the same `job_group` run as one cluster job, for example an annotation and the filters after it, per library:

```{graphviz}
digraph {
    rankdir=LR; bgcolor="transparent"; nodesep=0.2; ranksep=0.35;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b"];
    edge [color="#898781", arrowsize=0.6];
    calling [label="calling\nlibrary P001-T1"];
    subgraph cluster_group {
        label="one cluster job (job_group: tail)"; fontname="Helvetica"; fontsize=8; fontcolor="#52514e"; color="#e1e0d9"; style="rounded";
        annotation [label="annotation"]; filter1 [label="filtration_dkfz"]; filter2 [label="filtration_final"];
        annotation -> filter1 -> filter2;
    }
    calling -> annotation;
}
```

Their resource requests add up: runtimes of jobs that run one after another, cores and memory of jobs that run side by side.
So the members should share a partition.
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

`snappy report` builds the Snakemake report of the plots and QC outputs of the project's tasks.
It uses the same Snakefile, directory and profile as `snappy run`, and arguments after `--` go to Snakemake.
`--output` sets the file below the project directory.
A `.zip` report (the default `report.zip`) can hold HTML directories such as the FastQC reports.

### Contents

The report has one category per task (for example FastQC, Picard metrics, BAM statistics, `bcftools stats`, the copy number plots, peddy).
It embeds images, PDF and HTML; tables of statistics are downloads.
So a task whose tools write statistics also renders an HTML site into the report: a QC summary of all its libraries, or a calls report of its VCFs.

```{graphviz}
digraph {
    rankdir=LR; bgcolor="transparent"; nodesep=0.15; ranksep=0.4;
    node [shape=box, style="rounded,filled", fillcolor="#ffffff", color="#c3c2b7", fontname="Helvetica", fontsize=9, fontcolor="#0b0b0b"];
    edge [color="#898781", arrowsize=0.6];
    report [label="report.zip\none category per task", style="rounded,filled", fillcolor="#e1e0d9"];
    plots [label="plots and HTML reports\nembedded (PNG, PDF, HTML)"];
    tables [label="statistics tables\ndownloads"];
    qc [label="QC summary site\ntables and charts of all libraries"];
    calls [label="calls report site\nthe calls of each VCF, views, pedigrees"];
    igv [label="IGV views\nigv/top.html"];
    report -> plots; report -> tables; report -> qc; report -> calls; calls -> igv;
}
```

### QC summaries

A task whose tools write statistics renders a QC summary of all its libraries as tables and charts made with [datavzrd](https://datavzrd.github.io/).
The summary table is also a result file, `qc_summary/out/qc_summary.tsv`, with one row per library.

| Task                                            | Summary                                                                                                                                     |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| `ngs_mapping`                                   | samtools statistics: read counts and fractions per library, insert size and coverage distributions, base quality by cycle, reads per contig |
| `ngs_data_qc` with `tool: alfred`               | target coverage                                                                                                                             |
| `ngs_data_qc` with `tool: picard`               | one table per metrics class, the histograms as charts                                                                                       |
| `ngs_data_qc` with `tool: ngs_chew`             | the sample sheet checked against the fingerprints of all libraries, see below                                                               |
| `variant_calling` with `bcftools_stats` enabled | genotype counts, Ts/Tv, substitution, indel length, depth and quality distributions                                                         |
| `variant_checking`                              | peddy's checks of all pedigrees, see below                                                                                                  |

The site of a task opens with the summary table, whose cells show each value as a tick on the range of its column.
Its charts:

- The `libraries` chart puts one tick per library on each key metric, so a library apart from the others stands out and the tooltip names it.
- The distributions are thin lines, one per library, for up to six libraries (a color each, and the legend toggles them), and heatmaps of library rows for more.
- Charts zoom and pan, hovering a mark highlights its library, and a mark links to the library's row of the summary table.

The fingerprint summary of `ngs_data_qc` with `tool: ngs_chew` checks the sample sheet against the fingerprints of all libraries, one row per library:

- the sheet's sex against the sex predicted from the chrX genotypes (the sheet sex whose cohort median ratio is nearer, so the cohort needs libraries of both sexes),
- the sheet's father and mother against the relatedness and IBS0 rate of the pair,
- whether every pair of the pedigree has the relatedness the sheet implies (kinship over the pedigree, unsequenced donors included),
- the libraries of other pedigrees that are related at all.

Its charts are a scatter chart of the libraries by chrX ratio and chrY reads colored by the sheet's sex, a scatter chart of the pairs, and the table of the pairs.

The summary of `variant_checking` holds peddy's checks of all pedigrees: one row per sample with the sex check, heterozygosity, depth, call rate and ancestry, the pairs of every pedigree with their relatedness from the PED file and from the genotypes, and a scatter chart of the pairs on which a swapped, unrelated or duplicated sample falls off the diagonal.

### Calls report

A task that writes VCFs (`variant_calling`, `variant_filtration`, `variant_annotation`, `concat_variants`, `variant_phasing`, `combine_variants`) renders the calls of each VCF as a table when its `calls_report` section is set (`calls_report: {}` for the defaults).
The table is also the result file `report/calls/<entity>.calls.tsv`.

The columns follow what the VCF carries:

- Varlociraptor calls: the posterior probability of every event as a heatmap and per sample the allele frequency, the depth and the genotype.
  The details of a row show the observations (bars per strength of evidence, reads for the alternative allele to the right and for the reference to the left) and the posterior density of the allele frequency.
- Other callers: the genotype, allele frequency, alt reads, depth and genotype quality of each sample.
- Annotated variants (mehari's `ANN` or VEP's `CSQ`): the gene, consequence, impact and HGVS of the picked annotation; the frequency and ClinVar fields of mehari add the gnomAD frequency and the ClinVar classification.
- `linkouts` (gnomAD, ClinVar, dbSNP, Ensembl, GeneCards, OMIM and varlociraptor-inspect, which shows the observations of a call without any names) open in a new window.
  The ones that depend on the genome build need an annotation task in the chain, which tags its variants with its assembly.

The site has several views of the same calls:

- For varlociraptor calls, one view per event of the scenario (`denovo_<member>` for each child, `recessive`, `dominant`, `present`, or the somatic events) with the calls whose most probable event it is, ranked by that posterior, then gnomAD frequency (unobserved first), then impact.
  The site opens with the first view that has calls.
- Annotated calls of other callers get a `ranked` view by impact, frequency and quality.
- `all` lists every call in genomic order, without the per-sample plots and the inspect link, which the other views and the result TSV keep.
- A row links to the same call in the other views.

For a pedigree, the row details also draw the pedigree of the call:

- squares in aqua for male and circles in orange for female members (a legend under the pedigree), filled in that color when affected and white when not,
- dashed for unsequenced donors and unknown parents,
- a ring around the member of a per-member event such as `denovo_<member>`,
- for each member the call's allele frequency as a bar with its 90 % credible interval, its reads as a diverging bar (alternative allele to the right in blue, reference to the left in gray, darker for stronger evidence, on a gray background as wide as the most reads of a member), the depth and the genotype,
- the event with its posterior as the title.

With `calls_report: {igv: {}}` the task also reads the BAMs of `depends_on.alignments` and the genome and renders IGV views of the top calls of each view into one page, `igv/top.html` of the site.
`igv.top` (default 20) sets how many, round-robin over the views, `igv.flanking` the base pairs around a call, and `igv.max_size_mb` cuts the page.
The `igv` column of the table links a selected call to the page and every call to a running IGV desktop (`localhost:60151`).

## Log archive

When a run ends, successfully or not, `snappy run` writes `logs/<run start>.tar.gz`.
It holds:

- the Snakemake log of the run,
- the logs of the jobs that ran, with their SLURM logs,
- a copy of the `snkmt` database,
- in a git repository, `.snakemake/log/<run start>.git.txt`, which the start of the run writes with the commit of the configuration and the uncommitted changes of the tracked files with their diff (uncommitted changes are warned about).

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

The one exception is the BAM of an `ngs_mapping` task that a later task processes further, such as BQSR.
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

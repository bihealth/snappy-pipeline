(usage)=

# Usage

The `snappy` command is the interface to the pipeline.
`snappy --help` lists the subcommands.

(usage_commands)=

## Subcommands

The table lists each command with the first sentence of its help.
`snappy <command> --help` prints the full help.

```{snappy-commands}
```

The listing commands `snappy task list`, `snappy check-files` and `snappy rules` share one output style.
On a terminal they print a table; on a pipe, or with `--plain`, tab-separated lines without header and frame; with `--json`, an array with one object per row.
The data goes to stdout and the messages to stderr, so `snappy check-files --json > missing.json` holds only the rows.
`snappy scan-files` writes the files tables of the `external_data` tasks from the folders of the libraries (see {ref}`overview`); `snappy run --scan` runs it first.
`snappy sheet check` checks sample sheets and lists the cells that are synonyms of a value (`F` for `female`, `yes` for `true`, `exome` for `WES`), which `--fix` rewrites.
`snappy rules` lists the rules of the project as Snakemake names them: `<task>_<step>_<part>_<action>` for the jobs of a step part and `<task>_results` for the job that links the results of a task.
`--task` and `--part` restrict the listing, and `--names` prints only the names, one per line, so that a rerun of the QC summaries of every task is one command:

```shell
snappy run -- --forcerun $(snappy rules --part qc_summary --names)
```

`snappy task add` takes `[task_name=]step_type` arguments, for example `mapping=ngs_mapping`, and so does `snappy init --task`.
Without a task name, the task is named after its step.
On a terminal, `snappy task add` without arguments asks for the tasks one by one, and so does `snappy init` when it asks whether to add tasks: the step (type a part of its name to filter the list), the name, the tool and, for each `depends_on` key the tool requires, the task that provides it among those of the project, newest first.
For a key that no task provides yet, the question offers to add such a task on the spot, which pulls in a chain top-down: a calling task asks for its alignments, which asks for a mapping task, which asks for the genome, where an `external_data` task gets its `produces` and its file keys from the requirement.
The written block is that of the chosen tool: its section and the `depends_on` keys it reads, with the answers set and the keys still to fill in tagged `REQUIRED`, such as the paths of an `external_data` task.
Besides `config.yaml`, `snappy init`:

- puts the project into the directory given, else into the current directory if that is named like the project, else into a new directory named after the project.
- writes the header of an empty sample sheet `samplesheet.tsv`, which `config.yaml` names: the core columns and the ones a germline or cancer project usually fills (`family_id`, `father`, `mother`, `sex`, `is_affected`, `is_index`, `sample_name`, `is_tumor`, `library_kit`, `folder_name`); delete the columns the project does not need, add others (see the sample sheets section of the overview).
- records the project name, external id and description.
- makes the project a git repository that tracks only `.gitignore`, `config.yaml` and `README.md`.
  `--track-sheets` adds the sample sheet and the files table `reads.tsv`, `--remote <url>` adds the remote `origin` without pushing, and `--no-git` skips the repository.
- asks on a terminal for what no option gives.
  `--no-interactive` takes the defaults instead.

`--config-mode` sets the style of the generated task configuration.
`minimal`, the default, writes only the required fields, without descriptions.
`commented` writes every field with its description and comments out the optional ones, and `full` writes every field active.
`config.yaml` lists `pathvars`, `sample_sheets` and `resources` before `tasks`, so that added tasks end the file.
With `minimal`, the header of `config.yaml` points to the docs instead of explaining the sections.

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

On a terminal, `snappy run` offers to commit uncommitted changes to the tracked files of the project, such as `config.yaml`, and asks for the commit message.
Declined, without a terminal or in a dry-run, the run starts anyway: it records the commit and the diff and warns (see {ref}`Log archive <usage_log_archive>`).

Everything after `--` goes to Snakemake:

```shell
snappy run -- --dry-run --cores 1               # list the jobs, run nothing
snappy run --task mapping -- --cores 8          # one task and its missing upstream files
snappy run --task calling --frozen -- --cores 8 # reuse the existing mapping, never rerun it
snappy run -- --profile slurm                   # all tasks through the executor of a profile
```

## Cluster execution

`snappy run` always uses the workflow profile of snappy.
This profile deploys the software with conda, and with apptainer for the tools that run from a container image.
It also sets two loggers, `rich` for the console (styled output, progress bars with the pending jobs and the estimated remaining time per rule; `snappy run -- --logger snkmt` restores Snakemake's plain output) and `snkmt`, which `snappy watch` reads, and keeps Snakemake's metadata in a database instead of a directory of files.

The executor and its settings come from a Snakemake profile of the system or the user:

- below `/etc/xdg/snakemake/` for the system, or `~/.config/snakemake/` for the user,
- selected with `--profile <name>` after `--`, or with the environment variable `SNAKEMAKE_PROFILE`.

A SLURM profile, for example `~/.config/snakemake/slurm/config.yaml`:

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

On a cluster, every job runs as its own cluster job, unless its rule is local.
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

The resource requests of the grouped jobs add up: the runtimes of jobs that run one after another, the cores and memory of jobs that run side by side.
The jobs of a group should therefore share a partition.
`--group-components <group>=<n>` puts n independent groups into one cluster job, such as many small jobs per library:

```shell
snappy run -- --group-components tail=10
```

## Containers

PureCN and pVACtools run from container images, which their tool sections pin by digest (`purecn.container`, `pvactools.container`).
These fields also take the path of a local `.sif` file.
Snakemake pulls each image once into `.snakemake/software/` of the project and runs the jobs of these tools inside it with apptainer.
A job in a container sees `$HOME`, `/tmp`, the project directory and the snappy installation.
Add other directories that such a job reads, such as references outside the project, to `APPTAINER_BIND`, a comma-separated list:

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
: Holds only the result files.
    No task reads from it.
    It holds relative symlinks into `work/` with the same paths, and `output/.manifest.tsv`.
    For each task with result files, one local job, `<task>_results`, creates the links and the manifest.
    The manifest lists the path, size and md5 file of each result file.
    The job removes links to files that are no longer result files.
    A result file that is missing in `work/`, for example after a failed rerun, is built and linked again by the next `snappy run`.

Only result files have a `.md5` file next to them.
Logs are not result files.

## Report

`snappy report` builds the Snakemake report of the plots and QC outputs of the project's tasks.
It uses the same Snakefile, directory and profile as `snappy run`, and arguments after `--` go to Snakemake.
`--output` names the report file, relative to the project directory.
A `.zip` report (the default `report.zip`) can hold HTML directories such as the FastQC reports.

### Contents

The report has one category per task, with outputs such as FastQC, Picard metrics, BAM statistics, `bcftools stats`, the copy number plots and peddy.
It embeds images, PDF and HTML, and offers tables of statistics only as downloads.
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

The QC summary of a task shows the statistics of all its libraries as tables and charts made with [datavzrd](https://datavzrd.github.io/).
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

- The `libraries` chart puts one tick per library on each key metric.
  A library apart from the others stands out, and its tooltip names it.
- For up to six libraries, the distributions are thin lines, one per library in its own color.
  The legend toggles the lines.
  For more libraries, they are heatmaps with one row per library.
- Charts zoom and pan, hovering a mark highlights its library, and a mark links to the library's row of the summary table.

The fingerprint summary of `ngs_data_qc` with `tool: ngs_chew` checks the sample sheet against the fingerprints of all libraries, one row per library:

- the sheet's sex against the sex predicted from the ratio of heterozygous to homozygous chrX genotypes.
  A library gets the sheet sex whose cohort median ratio is nearer, so the cohort needs libraries of both sexes.
- the sheet's father and mother against the relatedness and IBS0 rate of the pair.
- whether every pair of a pedigree has the relatedness that the sheet implies, computed as kinship over the pedigree, unsequenced donors included.
- pairs of libraries from different pedigrees that are related at all.

The fingerprint summary has a scatter chart of the libraries by chrX ratio and chrY reads, colored by the sheet's sex.
It also has a scatter chart of the pairs and the table of the pairs.

The summary of `variant_checking` holds peddy's checks of all pedigrees.
It has one row per sample with the sex check, heterozygosity, depth, call rate and ancestry.
It lists the pairs of every pedigree with their relatedness from the PED file and from the genotypes.
A scatter chart shows the pairs, and a swapped, unrelated or duplicated sample falls off its diagonal.

### Calls report

A `variant_report` task renders the calls of each VCF of the task named in its `depends_on.variants` as a table.
Its `tool` is `datavzrd`, and any task that produces variants can be reported: calling, filtration, annotation, phasing, concatenated or combined calls.
The task writes one report per entity of the upstream task: per library, or per pedigree or donor when that task writes one VCF per related group.
The table is also the result file `report/calls/<entity>.calls.tsv`.

The columns follow what the VCF carries:

- Varlociraptor calls: the posterior probability of every event as a heatmap and per sample the allele frequency, the depth and the genotype.
  The details of a row show the observations and the posterior density of the allele frequency.
  The observations are bars per strength of evidence, with reads for the alternative allele to the right and for the reference to the left.
- Other callers: the genotype, allele frequency, alt reads, depth and genotype quality of each sample.
- Structural variants: `SV length` after the type, from `SVLEN` or else the span to `END`; a breakend has none.
- Annotated variants (mehari's `ANN` or VEP's `CSQ`): the gene, consequence, impact and HGVS of the picked annotation.
  The consequences are pills on one line, split at `&` in the annotator's order (mehari: the most severe first), each colored by its own impact and labelled without the `_variant` suffix, with the term on hover; `impact` is the highest of them.
  A call with more than three consequences shows two and a `+N` pill, which lists the rest on hover and shows them on click.
  The frequency and ClinVar fields of mehari add the gnomAD frequency and the ClinVar classification.
  Its CADD fields add `CADD` (the PHRED-scaled score, the raw score on demand), its SpliceAI fields `SpliceAI` (the largest delta score, the predictions per gene on demand), and its dbSNP annotation fills `id` with the rsID.
- `columns` sets the display mode of any column over these defaults, by its name in the result TSV: `normal`, `pinned`, `available`, `detail` or `hidden`, for example `columns: {CADD: pinned, gnomAD hom: normal}`.
  A column the table lacks is noted in the log.
- `linkouts` open in a new window: gnomAD, ClinVar, dbSNP, Ensembl, GeneCards, OMIM, the Sequence Ontology (each consequence pill) and varlociraptor-inspect.
  varlociraptor-inspect shows the observations of a call and puts no sample names into its link.
  The linkouts that depend on the genome build need an annotation task in the chain, because that task tags its variants with the assembly.
  Without such a task, they are left out.

The site has several views of the same calls:

- For varlociraptor calls, one view per event of the scenario: `denovo_<member>` for each child, `recessive`, `dominant`, `present`, or the somatic events.
  A view holds the calls whose most probable event it is, ranked by that posterior, then by gnomAD frequency (unobserved first), then by impact.
  The site opens with the first view that has calls.
- Annotated calls of other callers get a `ranked` view by impact, frequency and quality.
- `all` lists every call in genomic order, without the per-sample plots and the inspect link, which the other views and the result TSV keep.
- A row links to the same call in the other views.

For a pedigree, the row details also draw the pedigree of the call:

- Male members are squares in aqua and female members circles in orange, with a legend under the pedigree.
  A member is filled in its color when affected and white when not.
- Unsequenced donors and unknown parents are dashed.
- A ring marks the member of a per-member event such as `denovo_<member>`.
- Each member shows the call's allele frequency as a bar with its 90 % credible interval, its reads, the depth and the genotype.
  The reads are a diverging bar: the alternative allele to the right in blue, the reference to the left in gray, darker for stronger evidence.
  The gray background of this bar is as wide as the most reads of a member.
- The title is the event with its posterior.

With `igv: {}`, the task also renders IGV views of the top calls of each view into one page, `igv/top.html` of the site.
For these views, it reads the BAMs of `depends_on.alignments` and the genome of `depends_on.reference`.
`igv.top` (default 20) sets how many calls of each view get an IGV view, taken round-robin over the views.
`igv.flanking` sets the width of the window around a call in base pairs.
A structural variant spans to its `END`; one longer than `igv.max_span` (default 10 000 bp) gets a view of each breakpoint instead, so one large deletion cannot fill the page.
`igv.max_size_mb` limits the size of the page, which drops the lowest-ranked calls until the page fits.
The `igv` column of the table links a selected call to the page and every call to a running IGV desktop (`localhost:60151`).

(usage_log_archive)=

## Log archive

When a run ends, successfully or not, `snappy run` writes `logs/<run start>.tar.gz`.
It holds:

- the Snakemake log of the run,
- the logs of the jobs that ran, with their SLURM logs,
- a copy of the `snkmt` database,
- in a git repository, the record `.snakemake/log/<run start>.git.txt`.
  The run writes it at its start, with the commit of the configuration and the diff of uncommitted changes to tracked files.
  snappy warns about uncommitted changes.

Next to each job log lie the scripts of the job with every placeholder filled in: `<name>.sh`, the shell script as run, and `<name>.R` or `<name>.py` when the wrapper embeds an R or Python script.
From the project directory, `bash <name>.sh` reruns the job and `Rscript <name>.R` its R part; the archive holds the scripts.
snappy finds the job logs through the Snakemake log, so it does not search the task directories.
Dry-runs write no archive.

A killed Snakemake process (walltime, `scancel`) does not run its end hooks, so it writes no archive.
For such a run, `snappy logs` builds the archive from the newest Snakemake log, or from the log given as argument.

## Temporary files

Scatter fragments are the per-region or per-chunk pieces that a later rule of the task joins.
They are `temp()` files, and Snakemake deletes them when the joining job is done.

Per-sample intermediates are files that only later rules of the same task read, such as raw calls before filtering, coverage files, the normals of a panel and gVCFs.
They are kept, unless the task sets `delete_intermediates: [intermediates]`.

These files are never deleted:

- Files in `output/`.
- Files that a task passes to other tasks.
- Project-wide resources, such as installed software and prepared references.

So `--task` and `--frozen` runs find what they need.

The one exception is the BAM of an `ngs_mapping` or `alignment_processing` task that a later task processes further, such as BQSR.
With `delete_intermediates: [bam]`, the BAM and its index are temporary, have no `.md5` and are not linked into `output/`.
The BAM QC reports stay result files.
Snakemake deletes the BAM once the tasks that read it are done, and maps again when a later run needs it.
A `--frozen` run of a later task cannot map again, so it needs the BAM to still exist.

Once a temporary file is gone, a change to the rule that reads it (a parameter or its code) also reruns the rule that produced it.
For example, a change to the rule that joins the scatter fragments of a caller reruns all its scatter jobs.
With `delete_intermediates: [intermediates]`, a change to the filter also reruns the calling.

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

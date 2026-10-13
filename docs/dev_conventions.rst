.. _dev_conventions:

=====================
Developer Conventions
=====================

These rules apply to every workflow step. Tests enforce most of them; the section on each rule
names the check. ``tests/snappy_pipeline/test_conventions.py`` holds the static checks.

-------------------------
Configuration is explicit
-------------------------

A configuration must mean the same thing across snappy versions and setups. Code does not guess:
there is no automatic wiring of ``depends_on``, no fallback that picks a task by its step type, and
validation fails with a clear message instead of choosing for the user.

- ``tool`` has no default; every task names its tool (``test_tool_is_required``).
- ``depends_on`` keys have no task-name defaults (see `Contracts`_).

-------------
Step part API
-------------

The orchestrator loads the project once (``snappy_pipeline.orchestration.load_project``) and
creates the workflow object of every task once. A step Snakefile fetches it instead of building
its own:

.. code-block:: python

    from snappy_pipeline.orchestration import task_instance

    wf = task_instance(config["__task_name__"])

``wf.get_task_config(name)`` returns the task's own config (``name`` is its step or task name) or
the config of the upstream task set in ``depends_on.<name>``; anything else is an error.

Snakefiles stay thin and call the workflow object for everything. The methods a step part
provides:

.. list-table::
    :header-rows: 1

    * - Method
      - Returns
      - Rule section
    * - ``get_input_files(action)``
      - a function of wildcards, returning a dict (named inputs), a list or a path
      - ``input: unpack(wf.get_input_files(...))`` for dicts, otherwise ``wf.get_input_files(...)``
    * - ``get_params(action)``
      - a function of wildcards (optionally also of ``input``), returning a dict
      - ``params: args=wf.get_params(...)``
    * - ``get_output_files(action)``
      - static paths (Snakemake does not accept functions here)
      - ``output: **wf.get_output_files(...)``
    * - ``get_log_file(action)``
      - static paths
      - ``log: **wf.get_log_file(...)``
    * - ``get_resource(action, name)``
      - a function of wildcards, input, threads and attempt
      - ``threads:`` and ``resources:``

``BaseStep`` raises a ``TypeError`` at parse time if ``get_input_files`` or ``get_params`` returns
anything but a function.

Wrappers read the params dict as ``snakemake.params.args``. Snakemake does not allow
``unpack()`` in ``params:``, so every rule passes it as the single key ``args``.

Wrappers never read ``snakemake.config``. Every file a wrapper reads, such as the reference
FASTA, is a rule input, and every other value is a param. Only inputs order the jobs that
produce them first and trigger reruns when they change.
``test_wrappers_get_the_params_they_read`` checks that the jobs in the DAG snapshots pass every
``args`` key their wrapper reads.
``test_files_are_inputs_not_params`` flags params that hold file paths, against a baseline that
only shrinks.

Rules that run a wrapper declare ``benchmark: wf.get_benchmark_file(part, action)``: the path of
the rule's log below ``benchmarks/`` (``tasks/<task>/benchmarks/`` in a project), where Snakemake
writes the wall time, memory and I/O of each job. These files calibrate the resource requests.
Rules cached between workflows (``cache:``) cannot have a benchmark.
``test_wrapper_rules_write_a_benchmark`` and ``test_logs_and_benchmarks_are_unique_per_job``
check this.

A rule's log holds only its log file (plus a copy of the script, if the wrapper writes one).
Snakemake's metadata in ``.snakemake/`` records the code, params, inputs and conda environment
of every output, and the pinned environments (``environment.linux-64.pin.txt``) fix the package
versions.
Logs are not published; ``snappy run`` archives them per run.
Only published files, those linked into ``output/``, have an md5 file, and every published
file has one: SODAR needs one per uploaded file. A step declares ``<file>.md5`` next to
``<file>`` in the outputs, and the wrapper base classes write it. Wrappers never read md5 files.
A step's ``get_result_files`` lists the published files below ``output/``, each with its md5
file, and bgzipped VCF and BED files also with ``.tbi`` and ``.tbi.md5``. They mirror files
below ``work/``. One local job per task, ``<task>_publish`` (``snappy_pipeline/publish.py``),
links them and writes ``output/.manifest.tsv``; no other rule writes into ``output/``.
``get_output_paths`` returns ``work/`` paths, so tasks read each other's ``work/`` files and
never ``output/``.
``test_only_the_publish_job_writes_into_output``, ``test_no_job_reads_from_output``,
``test_logs_hold_only_logs``, ``test_logs_are_not_published``,
``test_md5_files_only_for_published_files``, ``test_published_files_come_with_md5_and_index``
and ``test_wrappers_read_no_md5_or_conda_files`` check this.

Every wrapper environment is pinned: Snakemake builds it from the exact package URLs in
``environment.linux-64.pin.txt`` next to ``environment.yaml``. After changing an environment
file, pin it again with ``uvx snakedeploy pin-conda-envs --conda-frontend conda <file>`` (this
needs conda). A symlinked ``environment.yaml`` gets a symlink to its target's pin file. The
``conda-autobump`` GitHub workflow updates and re-pins the environments weekly, with one pull
request per environment. ``test_wrapper_envs_are_pinned`` checks that each environment has its
pin file and that the pin file holds the exact versions of the environment file.

A step part lists the output keys of scatter fragments per action in ``temp_outputs``, and those
of per-sample intermediates that only later rules of the same task read in
``intermediate_outputs``. ``BaseStep.get_output_files`` marks fragments ``temp()`` always and
intermediates when the task sets ``delete_intermediates``. Published files, contract outputs,
files other tasks read and project-wide resources (installs, prepared references, pedigree
files) are neither. Every temp file needs a declared reader, otherwise Snakemake deletes it right
after its job; a job that reads a temp file's index implicitly declares the index as input. The
DAG snapshots are built with ``delete_intermediates: true``, and ``test_intermediates_are_temp``
and ``test_temp_files_are_intermediates`` check them.

Outputs for the Snakemake report are declared per action in a step part's ``report_outputs``
(output key to ``ReportOutput``: subcategory, caption, ``htmlindex`` for a directory).
``BaseStep.get_output_files`` flags them with ``report()``, in the task's category and labelled
by the wildcards of the path. Captions are RST files in ``report/`` next to the step module.

Defining input files and params
===============================

``BaseStepPart.get_input_files(action)`` and ``get_params(action)`` validate the action against
the class's ``actions`` tuple and return ``self._get_input_files_<action>`` or
``self._get_params_<action>``. A step part therefore only defines these methods:

.. code-block:: python

    class ExampleStepPart(BaseStepPart):
        name = "example"
        actions = ("run",)

        def _get_input_files_run(self, wildcards):
            return {"vcf": "work/{library_name}/out/{library_name}.vcf.gz".format(**wildcards)}

        def _get_params_run(self, wildcards):
            return {"min_af": self.config.example.min_af}

Rules:

- ``actions`` is a tuple.
- Input and params methods take exactly ``(self, wildcards)``; params may also take
  ``(self, wildcards, input)``. No ``**kwargs``: Snakemake passes every available argument to a
  function that accepts ``**kwargs``.
- Override ``get_input_files`` or ``get_params`` only if it does more than dispatching, and then
  still return a function.
- ``@dictify`` and ``@listify`` keep the decorated function's signature (``functools.wraps``),
  which Snakemake inspects.
- Paths returned by input functions are relative (``work/...``, or upstream paths from
  ``get_upstream_paths()``). Snakemake adds the task prefix to function results just as it does
  to static patterns.

The exception is an input function that needs Snakemake's ``checkpoints`` object; the Snakefile
calls it as ``wf.get_input_files(...)(wildcards, checkpoints)`` from a Snakefile-level function
(``helper_gcnv_model_*``).

---------
Contracts
---------

A task's outputs are described by ``DataSignature`` objects: a data type plus tags such as
``dna``, ``somatic`` or ``filtered``.

- A step lists what its tasks produce in ``produces``. If that depends on the task's config or
  its inputs, the step overrides the classmethod ``task_produces(config, upstream)``. For
  example, ``ngs_mapping`` produces ``alignments [rna]`` for STAR, and ``variant_filtration``
  adds ``filtered`` to the tags of its input.
- ``depends_on`` keys are named after the data they bring in, not after the step that produces
  it: ``reads``, ``alignments``, ``variants``, ``copy_number``, ``structural_variants``,
  ``fusions``, ``expression``, ``hla_types``, ``strandedness``, ``panel_of_normals``,
  ``reference``, ``features``, ``dbsnp`` and ``index``. A step with two inputs of one kind
  qualifies them by role, such as ``somatic_variants`` and ``germline_variants``.
  ``test_depends_on_keys_use_the_vocabulary`` checks this.
- ``reference`` (genome FASTA), ``features`` (gene annotation GTF) and ``dbsnp`` use the
  shared field types ``Reference``, ``Features`` and ``Dbsnp`` from
  ``snappy_pipeline.workflows.abstract.protocol``. They usually name ``external_data`` tasks.
- A key that every tool of a step reads is required and has no default. A key that only some
  tools read defaults to ``""`` and is listed per tool in the step's ``TOOL_DEPENDENCIES``,
  which ``validators.require_tool_dependencies`` checks. Other optional keys default to ``""``,
  meaning "not used". No key defaults to a task name.
- ``reads`` names an ``external_data`` task, a task that writes FASTQs (such as
  ``adapter_trimming``), or the reserved value ``data_sets``, which searches the FASTQs in the
  data sets' search paths. No task may be called ``data_sets``.
- Each ``depends_on`` field states what it requires with a ``DataSignature`` in its
  ``Annotated`` metadata. This annotation is the only place a requirement is declared.
- A task can produce several signatures. mutect2 calling produces all calls and the PASS
  calls, which carry the ``filtered`` tag. ``select_signature`` picks the one a consumer reads:
  the ``filtered`` one when several satisfy the requirement, so PASS calls are the default and
  ``-filtered`` asks for all calls. ``get_output_paths`` receives the picked signature and
  returns the matching file. ``get_upstream_paths(field, signature=...)`` narrows the
  requirement for a single lookup.
- ``load_project()`` computes the signatures of every task in dependency order and checks
  each requirement before Snakemake builds any rule. A mismatch fails with a message that names
  both tasks.
- Consumers read upstream files only through ``get_upstream_paths(field, library_name=...)``,
  which returns the provider's named outputs (``.bam``, ``.vcf``, ...) wrapped in the field's
  ``ExpectedPathSchema``; the schema receives only the keys it declares. Consumers never build a
  provider's paths themselves (``test_upstream_paths_come_from_contracts``).
- ``get_output_paths(config, signature, ...)`` receives the config of the task whose outputs
  are requested. It does not check the signature again; it only uses it to choose between
  several outputs.

-----
Paths
-----

- Only truly variable values are wildcards: libraries, samples, tumor and normal, cohorts.
- The tool is part of the task configuration and never appears in a path, neither as a directory
  prefix (``work/cnvkit.{library_name}/``) nor in a file name.
- Files are named ``{entity}.{suffix}``, with a dot after the wildcard, e.g.
  ``{tumor_library}.gene_log2.txt``. Underscores stay inside multi-word suffixes.

``test_path_conventions_only_improve`` checks the last two rules on all DAG snapshots against a
baseline of existing violations (``tests/snappy_pipeline/snapshots/path_convention_violations.json``).
New violations fail the test; fixed ones must be removed from the baseline.

--------------------
I/O footprint on HPC
--------------------

On shared file systems such as cephfs, and with slurm, the number of small files, metadata
operations (stat, readlink) and jobs costs more than bytes moved. Keep all three low: no extra
jobs for bookkeeping, no companion files without a consumer, no declared outputs nobody needs.

-------
Testing
-------

``tests/snappy_pipeline/test_generated_configs_dryrun.py`` builds the DAG of every generated task
closure and compares each job (rule, wildcards, input, output and log paths, params, threads,
resources and wrapper) with a snapshot in ``tests/snappy_pipeline/snapshots/dag/``:

- A refactoring must leave the snapshots unchanged.
- An intended change is reviewed as a snapshot diff; refresh the snapshots with
  ``SNAPPY_UPDATE_SNAPSHOTS=1`` (see the contributing guide).
- Every closure must target at least one file. Exceptions are listed with a reason in
  ``EXPECTED_EMPTY``; closures that cannot be built yet are strict xfails in ``KNOWN_BROKEN``.

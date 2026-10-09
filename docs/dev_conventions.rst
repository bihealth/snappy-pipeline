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

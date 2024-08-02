.. _dev_intro:

========================
Developer's Introduction
========================

.. note::

    Before reading this chapter, you should

    - have knowledge from the user's perspective of CUBI pipeline (start a :ref:`usage`).

    After reading this chapter, you should

    - know about the Python programming techniques required from a CUBI pipeline developer
    - have an overview of the components of a pipeline step
    - know that ``cubi-snake`` only serves as a shortcut to the ``snakemake`` executable.

The target audience of this part of the documentation is developers who want to change or extend the pipeline.
The aim is to give a good overview of the architecture of the pipeline system and dissect some typical existing pipeline steps for educational purposes.
Most parts of the system follow a consistent programming and architecture style that should be followed to ease the understanding of the system.

If you are a proficient Python programmer then you should not have a too hard time to get started.
If your Python karate is less strong (e.g., if you are a Bioinformatician coming from the "bio" and not the "informtician" side), take a deep breath and brace yourself, you will learn something here.
Before we start, here is the Zen of Python as a reminder::

    >>> import this
    The Zen of Python, by Tim Peters

    Beautiful is better than ugly.
    Explicit is better than implicit.
    Simple is better than complex.
    Complex is better than complicated.
    Flat is better than nested.
    Sparse is better than dense.
    Readability counts.
    Special cases aren't special enough to break the rules.
    Although practicality beats purity.
    Errors should never pass silently.
    Unless explicitly silenced.
    In the face of ambiguity, refuse the temptation to guess.
    There should be one-- and preferably only one --obvious way to do it.
    Although that way may not be obvious at first unless you're Dutch.
    Now is better than never.
    Although never is often better than *right* now.
    If the implementation is hard to explain, it's a bad idea.
    If the implementation is easy to explain, it may be a good idea.
    Namespaces are one honking great idea -- let's do more of those!


-------------------------------
Prerequisites -- Your Tool Belt
-------------------------------

The CUBI pipeline system is implemented using Python 3 (>=3.4 at the moment) and built upon the wonderful `Snakemake <https://snakemake.bitbucket.org>`_ (>=3.10 at the moment).
For distributed, parallel execution, the pipeline is tailored towards execution with SGE Grid Engine.
In order to follow this developer's documentation comfortably, you should be familiar with all three systems:

- Python 3
- Snakemake
- Grid Engine (or similar cluster job queueing system).

You should be familiar with the CUBI pipeline from the user perspective already.

Also, understanding in the following techniques will come in handy:

- Snakemake as a rule-based language for describing workflows,
- object oriented programming,
- Python generators (via ``yield``, use of ``yield from``),
- Python decorators,
- Python ``itertools`` package and built-ins such as ``zip``/``map``
- Exception handlng
- JSON, JSON schema,
- scope, lambdas, and closures in Python,
- realize that classes are objects themselves and callable (their constructor),
- the ``snakemake.io.expand()`` helper function
- text manipulation using ``str.format``, ``textwrap.dedent``, ``str.lstrip``,
- understanding in common Python standard library code such as ``os[.path]``, ``collections.OrderedDict``,
- the recently added Snakemake ``unpack()`` keyword,
- the concept of mixin classes.

The following are handy references about using Python effectively:

- `The Hitchhiker's Guide to Python <http://docs.python-guide.org/en/latest/>`_.
- Slatkin, Brett: Effective Python: 59 Specific Ways to Write Better Python


----------------------------------
Anatomy of a Typical Pipeline Step
----------------------------------

Each pipeline step is implemented as a Snakemake workflow.
For each step, there is a module sub directory below ``snappy_pipeline.workflows`` containing:

- ``__init__.py`` with classes that actually implement the workflow
- ``Snakefile`` that contains the Snakemake rule definitions but usually just hooks in calls to the actual implementation code from ``__init__.py``.

Usually, you define a :class:`BaseStep <snappy_pipeline.workflows.abstract.BaseStep>` sub class in your Python code (``__init__.py``) that is then instantiated in your ``Snakefile``.
The current configuration is passed into the constructor of this class and it then "takes over" and applies default setting, generating cluster resource settings, etc.
Then, you pass the result of method calls to your :class:`BaseStep <snappy_pipeline.workflows.abstract.BaseStep>` instance as the values for the ``input:``, ``output:``, etc. sections of your ``Snakefile``.

.. warning::
   By convention your new Workflow step should be instantiated as ``wf = StepClass(...)`` in the ``Snakefile`` during object setup. Otherwise tools including cubi-tk might not be able to detect and parse your step. See existing workflow ``Snakefile`` for reference.

The :class:`BaseStep <snappy_pipeline.workflows.abstract.BaseStep>` sub class itself uses :class:`BaseStepPart <snappy_pipeline.workflows.abstract.BaseStepPart>` sub classes for the implementation of the individual parts.
One part might be linking in FASTQ files from the raw input directory or linking from the ``work/`` to the ``output/`` directory.
Another part might be the somatic variant calling using mutect or WGS SV calling using Delly2.

Each of the parts might be split into different actions if the implementing tools need their own more or less complex "workflow" themselves.
An example for such a tool is Delly2 where first variant calling is performed for each sample, then the resulting site list is merged and used for genotpying is all samples individually.
Finally, the wohle cohort's genotypes are merged and for each sample, only the variants that have been observed in it will be executed.
If the tools can just be executed in one action, this action should be called ``"run"``.

This approach has the advantage that most complex things happen in Python code for which tools for testing, (some) static code analysis, documentation, and style checking exist.
In the Python files, we can use the whole Python tooling ecosystem whereas in the Snakemake files, tools would choke on the first ``rule`` keyword.
In short, the ``Snakefile`` only serves as the entry point for your Python code.


----------------------------------------
Anatomy of the ``cubi-snake`` Executable
----------------------------------------

CUBI pipeline runs are invoked with the ``cubi-snake`` executable that internally calls Snakemake with sensible defaults for either local execution or execution on via SGE on an HPC cluster.
It serves as a convenience wrapper that reads the current pipeline step from the current working directories ``config.yaml`` file (where available, otherwise you have to use the ``--step`` argument).

Some parameters are handed through directly to Snakemake, others are serve as macros that add more complex parameters with best pratice values or print the configuration setting.

This sounds like an aweful amount of "magic" but is quite simple and transparent, really.
The generally useful ``snakemake`` parameters are also available to ``cubi-snake`` (or should be added, please create a ticket).
Also, snakemake is invoked through the command line interface and a command line to copy and paste is printed at the beginning of every ``cubi-snake`` invocation.

.. _contributing:

.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every
little bit helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/bihealth/snappy-pipeline/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug"
and "help wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

CUBI Pipeline could always use more documentation, whether as part of the
official CUBI Pipeline docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/bihealth/snappy-pipeline/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `snappy-pipeline` for local development.

1. Fork the `snappy_pipeline` repo on BIH GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:bihealth/snappy-pipeline.git

3. Install the development environment with `pixi <https://pixi.sh>`_ (Python 3.12 or newer)
   and the pre-commit hooks (ruff and snakefmt)::

    $ cd snappy-pipeline/
    $ pixi install --environment dev
    $ pixi run -e dev pre-commit install

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, format the code and check that the linters and the tests
   pass::

    $ pixi run -e dev srcfmt
    $ pixi run -e dev lint
    $ pixi run -e dev test

6. Commit your changes and push your branch to GitHub. Use
   `Conventional Commits <https://www.conventionalcommits.org>`_ for commit messages and pull
   request titles (e.g. ``fix(cli): ...``, ``feat(ngs_mapping): ...``); CI checks the pull
   request title::

    $ git add <changed files>
    $ git commit -m "fix(scope): short description of the change"
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring.
3. Linting and tests pass in CI (``.github/workflows/ci.yml``). Changes to workflows or
   wrappers are also dry-run and run end to end on the pipelines in ``.tests/``
   (``.github/workflows/ci-e2e.yml``).

Tips
----

To run a subset of tests::

$ pixi run -e dev test tests/snappy_pipeline/apps

To dry-run one of the test pipelines::

$ cd .tests/test-workflow/pipelines/snappy-germline_wes
$ pixi run snappy run -n -- --cores 1

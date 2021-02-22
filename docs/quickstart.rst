.. _quickstart:

==========
Quickstart
==========

This chapter gives the minimal number of commands required for setting you up with the pipeline on the BIH cluster with DRMAA.

.. note::

    This describes the setup as a pipeline user.
    If you want to know about the setup as a pipeline developer, see :ref:`installation`.

-------------------
Install (Mini)conda
-------------------

First, install miniconda, e.g., into ``$HOME/miniconda3``.

.. code-block:: shell

    $ wget -O /tmp/Miniconda3-latest-Linux-x86_64.sh \
        https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ bash /tmp/Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3

.. note:: What is conda/miniconda?

    Conda is a Python-based package manager that can also package binary files (such as Bioinformatics software).
    Miniconda is a minimal Conda installation.

    If anything goes wrong with your Miniconda installation, you can always just remove ``$HOME/miniconda3`` and start anew.

Now, make sure it is available in your ``PATH`` environment variable.

.. code-block:: shell

    $ export PATH=$HOME/miniconda3/bin:$PATH

-----------
DRMAA Setup
-----------

Setup the environment variable in your ``~/.bashrc`` file and **then** reconnect (open a new connection) to a cluster node with ``qrsh``.
It is important to have ``DRMAA_LIBRARY_PATH`` defined before installation as this will install the Python ``drmaa`` package which is necessary for robust job submission to the cluster.

You can skip this step if you have everything setup properly in your ``~/.bashrc``.

This command will extend your ``~/.bashrc`` file to properly use your Miniconda setup.

.. code-block:: shell

    $ cat <<"EOF" >>~/.bashrc
    case "${HOSTNAME}" in
        med-login*)
            ;;
        med-transfer*)
            ;;
        *)
            # properly setup TMPDIR
            export TMPDIR=$HOME/scratch/tmp
            mkdir -p ${TMPDIR}

            # Make Miniconda3 installation available
            export PATH=$HOME/miniconda3/bin:$PATH

            # export DRMAA library
            export DRMAA_LIBRARY_PATH=/opt/sge/lib/lx-amd64/libdrmaa.so
            ;;
    esac
    EOF

Now, make sure to reconnect to a cluster node.
The ``DRMAA_LIBRARY_PATH`` environment variable must be set:

.. code-block:: shell

    $ qrsh
    [...]
    $ echo $DRMAA_LIBRARY_PATH
    /opt/sge/lib/lx-amd64/libdrmaa.so

---------------------
Install CUBI Pipeline
---------------------

The recommended way of installing the CUBI pipeline is via ``pip``.

Replace the `X.Y.Z` in the definition of ``VERSION`` below with the version you find in the ``README.rst`` file of the project on the CUBI Gitlab.

.. code-block:: shell

    $ VERSION=vX.Y.Z
    $ pip install git+ssh://git@gitlab.bihealth.org/cubi/snappy_pipeline.git@v${VERSION}#egg=snappy_pipeline

Execute ``cubi-snake --cubi-pipeline-self-test`` to see the current setup.
Make sure that the configuration shows that DRMAA is enabled.

.. code-block:: shell
    :emphasize-lines: 10

    $ cubi-snake --cubi-pipeline-self-test
    CUBI Pipeline
    =============

    Version     <installed_version>

    Features
    --------

    DRMAA       yes

.. note:: What is DRMAA?

    **DRMAA** (DRMAA or Distributed Resource Management Application API) is an interface for the **submission** and **control** of jobs to cluster schedulers, such as Sun Grid Engine (SGE) or Univa Grid Engine (UGE).
    Submission with DRMAA is similar to using the ``qsub`` command, but using it allows Snakemake (and thus the CUBI pipeline) to track the submitted jobs and detect failures instead of waiting indefintely.

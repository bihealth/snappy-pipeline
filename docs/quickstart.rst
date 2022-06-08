.. _quickstart:

==========
Quickstart
==========

This chapter gives the minimal number of commands required for setting up the pipeline on the BIH cluster.

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

---------------------
Install CUBI Pipeline
---------------------

The recommended way of installing the CUBI pipeline is via ``pip``.

Replace the `X.Y.Z` in the definition of ``VERSION`` below with the version you find in the ``README.rst`` file of the project on the CUBI GitHub.

.. code-block:: shell

    $ VERSION=vX.Y.Z
    $ pip install git+ssh://git@github.com:bihealth/snappy-pipeline.git@v${VERSION}#egg=snappy_pipeline

Or see ``README.rst`` for a more detailed installation guide and the environment setup step.

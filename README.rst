===============================================
SNAPPY - SNAPPY Nucleic Acid Procesing Pipeline
===============================================

- **Version**: ``v0.3.0``
- **License**: MIT

------------
Installation
------------

**In a nutshell**:

.. code-block:: bash

    # Download & preparation
    git clone git@github.com:bihealth/snappy-pipeline.git
    cd snappy-pipeline

    # If you want to select a given branch, uncomment the following:
    # git checkout <branch_name>

    # WARNING- make sure that you are in your conda base environment

    # Create conda environment "snappy_env" with the minimal requirements:
    mamba env create --file environment.yml
    conda activate snappy_env
    pip install --file requirements/drmaa.txt

    # Add testing & development requirements:
    pip install --file requirements/test.txt
    pip install --file requirements/dev.txt

    # Optionally add "pytest-pdb" missing from anaconda
    pip install pytest-pdb

    # Install snappy in snappy_env environment
    pip install -e .

Installation should be complete in 10 to 15 minutes.

**Note:** To create the environment under another name, replace the commands for the environment creation & activation of the correct environment by:

.. code-block:: bash

    mamba env create --file environment.yml --name <other_environment_name>
    conda activate <other_environment_name>


See `user installation <docs/quickstart.rst>`_ if you just want to use the pipeline.

See `developer installation <docs/installation.rst>`_ for getting started with working on the pipeline code and also building the documentation.


#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Installation driver (and development utility entry point) for snappy-pipeline
"""

from itertools import chain
import os
import sys

from setuptools import find_packages, setup


__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def parse_requirements(path):
    """Parse ``requirements.txt`` at ``path``."""
    requirements = []
    with open(path, "rt") as reqs_f:
        for line in reqs_f:
            line = line.strip()
            if line.startswith("-r"):
                fname = line.split()[1]
                inner_path = os.path.join(os.path.dirname(path), fname)
                requirements += parse_requirements(inner_path)
            elif line != "" and not line.startswith("#"):
                requirements.append(line)
    return requirements


# Enforce python version >=3.7
if sys.version_info < (3, 7):
    print("At least Python 3.7 is required.\n", file=sys.stderr)
    sys.exit(1)

with open("README.md") as readme_file:
    readme = readme_file.read()

with open("CHANGELOG.md") as history_file:
    history = history_file.read()

# Get requirements
requirements = parse_requirements("requirements/base.txt")

test_requirements = [
    # TODO: put package test requirements here
]


# Name of the tools
TOOLS = (
    "bed_filter_jaccard",
    "bam_cov_stats",
    "fix_vcf",
    "genome_windows",
    "quickvenn",
    "vcf_first_header",
    "vcf_filter_denovo",
    "vcf_filter_from_info",
    "vcf_filter_to_info",
    "ped_to_vcf_header",
)

# Name of the scripts
SCRIPTS = ("snappy-vcf_sort", "snappy-vcf_split")


def console_scripts_entry_points(names, module):
    """Yield entries for the 'console_scripts' entry points"""
    prefix = "snappy"
    for name in names:
        if module == "wrappers":
            yield "{prefix}-{name} = snappy_wrappers.{module}.{name}.__main__:main".format(
                prefix=prefix, name=name, module=module
            )
        elif module == "tools":
            yield "{prefix}-{name} = snappy_wrappers.{module}.{name}:main".format(
                prefix=prefix, name=name, module=module
            )


def bash_scripts(names):
    """Return generator with bash scripts of the given ``names``"""
    return (os.path.join("scripts", name) for name in names)


package_root = os.path.abspath(os.path.dirname(__file__))

version = {}
with open(os.path.join(package_root, "snappy_pipeline/_version.py")) as fp:
    exec(fp.read(), version)
version = version["__version__"]

setup(
    name="snappy-pipeline",
    version=version,
    description="SNAPPY Nucleic Acid Processing in Python (by CUBI)",
    long_description=readme + "\n\n" + history,
    long_description_content_type="text/markdown",
    author="Manuel Holtgrewe",
    author_email="manuel.holtgrewe@bih-charite.de",
    url="https://github.com/bihealth/snappy-pipeline",
    packages=find_packages(),
    package_dir={"snappy_wrappers": "snappy_wrappers", "snappy_pipeline": "snappy_pipeline"},
    entry_points={
        "console_scripts": list(
            chain(
                (
                    "snappy-pull-sheet = snappy_pipeline.apps.snappy_pull_sheet:main",
                    "snappy-snake = snappy_pipeline.apps.snappy_snake:main",
                    "snappy-start-project = snappy_pipeline.apps.snappy_start_project:main",
                    "snappy-start-step = snappy_pipeline.apps.snappy_start_step:main",
                    "snappy-refresh-step = snappy_pipeline.apps.snappy_refresh_step:main",
                    "snappy-slurm-status = snappy_pipeline.apps.snappy_slurm_status:main",
                ),
                console_scripts_entry_points(TOOLS, "tools"),
            )
        )
    },
    scripts=tuple(bash_scripts(SCRIPTS)),
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords="bioinformatics",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
        # We are missing bioconda pysam packages for 3.11 and 3.12, cf.
        # https://github.com/bioconda/bioconda-recipes/issues/37805
        # "Programming Language :: Python :: 3.11",
        # "Programming Language :: Python :: 3.12",
    ],
    test_suite="tests",
    tests_require=test_requirements,
)

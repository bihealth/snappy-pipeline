# -*- coding: utf-8 -*-
"""Code for running a sequential CUBI+Snakemake wrappers
"""

import json
import os
import pickle
import sys
import tempfile

from snakemake import shell, snakemake
from snakemake.script import Snakemake
from termcolor import colored

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


#: Template for the temporary Snakemake file
SNAKEMAKE_TPL = r"""
# Auto-generated CUBI+Snakemake Snakemake file for using Snakemake wrappers
import pickle; obj = pickle.loads({pickled_obj})

configfile: {path_config}

rule all:
    input: '.wrapper.done'

rule main:
    output: touch('.wrapper.done'){output_token}
    input: {prefixes[input]}obj.input
    params: {prefixes[params]}obj.params
    threads: obj.threads
    resources: {prefixes[resources]}obj.resources
    log: obj.log
    wrapper: obj.params['__wrapper']
""".lstrip()


def run(wrapper_info):
    """Run sequential CUBI+Snakemake wrapper from given :py:class:`~snappy_wrappers.base.WrapperInfo`
    object
    """
    shell.executable("/bin/bash")
    args = wrapper_info.args
    with tempfile.TemporaryDirectory() as tmpdirname:
        obj = Snakemake(
            wrapper_info.get_input(),
            wrapper_info.get_output(),
            wrapper_info.get_params(),
            None,
            wrapper_info.get_threads(),
            wrapper_info.get_resources(),
            wrapper_info.get_log(),
            wrapper_info.get_config(),
            str(wrapper_info.__class__),
        )
        path_config = os.path.join(tmpdirname, "config.json")
        with open(path_config, "wt") as fconfig:
            json.dump(wrapper_info.get_config(), fconfig)
        contents = SNAKEMAKE_TPL.format(
            pickled_obj=repr(pickle.dumps(obj)),
            prefixes=wrapper_info.get_value_prefixes(),
            path_config=repr(path_config),
            output_token=(
                ", {}obj.output".format(wrapper_info.get_value_prefixes()["output"])
                if obj.output
                else ""
            ),
        )
        path_snakefile = os.path.join(tmpdirname, "Snakefile")
        with open(path_snakefile, "wt") as f:
            f.write(contents)
        if args.verbose:
            print(
                colored("Creating temporary working directory", "yellow") + " " + tmpdirname,
                file=sys.stderr,
            )
            print(colored("Snakemake file contents", "yellow"), file=sys.stderr)
            print(contents, file=sys.stderr)
            print(colored("Launching Snakemake...", "blue"), file=sys.stderr)
        snakemake(
            path_snakefile,
            workdir=tmpdirname,
            verbose=args.verbose,
            quiet=not args.verbose,
            use_conda=True,
            overwrite_shellcmd="/bin/bash",
        )

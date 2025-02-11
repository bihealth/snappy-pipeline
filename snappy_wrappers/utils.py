# -*- coding: utf-8 -*-
"""Utility code for snappy_wrappers"""

import os
import re
import subprocess
import json

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"


def basename(fname, strip_suffix=""):
    """Unix `basename`-like function"""
    fname = os.path.basename(fname)
    if strip_suffix and fname.endswith(strip_suffix):
        return fname[: -len(strip_suffix)]
    else:
        return fname


def listify(gen):
    """Decorator that converts a generator into a function which returns a list

    Use it in the case where a generator is easier to write but you want
    to enforce returning a list::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield i
    """

    def patched(*args, **kwargs):
        """Wrapper function"""
        return list(gen(*args, **kwargs))

    return patched


def dictify(gen):
    """Decorator that converts a generator into a function which returns a dict

    Use it in the case where a generator is easier to write but you want
    to enforce returning a dict::

        @listify
        def counter(max_no):
            i = 0
            while i <= max_no:
                yield 'key{}'.format(i), i
    """

    def patched(*args, **kwargs):
        """Wrapper function"""
        return dict(gen(*args, **kwargs))

    return patched


def install_R_package(
    dest: str, name: str, repository: str = "cran", url: str | None = None) -> subprocess.CompletedProcess:
    assert dest, "Missing R package destination folder"
    os.makedirs(os.path.dirname(dest), mode=0o750, exist_ok=True)

    match repository:
        case "cran":
            install_cmd = f"install.packages('{name}', lib='{dest}', repos='https://cloud.r-project.org', update=FALSE, ask=FALSE)"
        case "bioconductor":
            install_cmd = f"BiocManager::install('{name}', lib='{dest}', update=FALSE, ask=FALSE)"
        case "github":
            assert url, f"Can't install R package '{name}' from github, URL is missing"
            install_cmd = f"remotes::install_github('{url}', lib='{dest}', upgrade='never')"
        case "bitbucket":
            
            print("#############################################################################################################################DEBUGGING#############################################################")
            assert url, f"Can't install R package '{name}' from bitbucket, URL is missing"
            install_cmd = f"remotes::install_bitbucket('{url}', lib='{dest}', upgrade='never')"
        case "local":
            assert url, f"Can't install local R package '{name}', missing path"
            assert os.path.exists(url), f"Can't find local R package '{name}' at location '{url}'"
            install_cmd = f"install.packages('{url}', repos=NULL, lib='{dest}', update=FALSE, ask=FALSE)"
        case _:
            raise ValueError("Unknown repository '{repository}'")
    R_script = [
        f".libPaths(c(.libPaths(), '{dest}'))",
        install_cmd,
        f"status <- try(find.package('{name}', lib.loc='{dest}', quiet=FALSE, verbose=TRUE))",
        "status <- ifelse(is(status, 'try-error'), 1, 0)",
        "quit(save='no', status=status, runLast=FALSE)",
    ]
    cmd = ["R", "--vanilla", "-e", "; ".join(R_script)]
    return subprocess.run(cmd, text=True, check=True)


def install_R_packages(dest: str, filename: str):
    with open(filename, "rt") as f:
        packages = json.load(f)
    for package in packages:
        status = install_R_package(dest, name=package["name"], repository=package["repository"], url=package.get("url", None))
        status.check_returncode()

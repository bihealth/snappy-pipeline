# -*- coding: utf-8 -*-
"""Utility code for snappy_wrappers"""

import os
import re
import subprocess

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


def install_R_package(dest: str, name: str, repo: str):
    """Installs R package <name> in directory <dest>

    <repo> can be on of "cran", "bioconductor", "github", "bitbucket" or "local".

    Github & bitbucket packages can be further defined in the <name>,
    choosing a commit, a tag or a pull request (see the remotes package reference).
    For local packages, the <name> must be the path to the package source.

    In these cases, the routine tries to be clever and guess the package name from
    the <name>. It implements the following recipes:

    - basename of the path without (tar/zip) extension for local packages
    - remove username, subdirectory and reference/release/pull request for github & bitbucket

    When installation isn't successful, a subprocess.CalledProcessError is raised
    """

    assert repo in (
        "cran",
        "bioconductor",
        "github",
        "bitbucket",
        "local",
    ), f"Unknown/unimplemented repository {repo}"

    os.makedirs(os.path.dirname(dest), mode=0o750, exist_ok=True)

    if repo == "cran":
        install_cmd = f"install.packages('{name}', lib='{dest}', update=FALSE, ask=FALSE)"
    elif repo == "bioconductor":
        install_cmd = f"BiocManager::install('{name}', lib='{dest}', update=FALSE, ask=FALSE)"
    elif repo == "github":
        path = name
        pattern = re.compile("^([^/]+)/([^/]+)(/[^@|]+)?([@|].+)?$")
        m = pattern.match(path)
        assert m, f"Cannot extract package name from github path {path}"
        name = m.groups()[1]
        install_cmd = f"remotes::install_github('{path}', lib='{dest}', upgrade='never')"
    elif repo == "bitbucket":
        path = name
        pattern = re.compile("^([^/]+)/([^/]+)([/@].+)?$")
        m = pattern.match(path)
        assert m, f"Cannot extract package name from bitbucket path {path}"
        name = m.groups()[1]
        install_cmd = f"remotes::install_bitbucket('{path}', lib='{dest}', upgrade='never')"
    elif repo == "local":
        path = name
        pattern = re.compile("^(.+?)(\\.(zip|tar(\\.(gz|bz2))?|tgz2?|tbz))?$")
        m = pattern.match(os.path.basename(path))
        assert m, f"Cannot extract package name from local filesystem path {path}"
        name = m.groups()[0]
        install_cmd = f"remotes::install_local('{path}', lib='{dest}', upgrade='never')"
    else:
        install_cmd = None

    R_script = [
        install_cmd,
        f"status <- try(find.package('{name}', lib.loc='{dest}', quiet=FALSE, verbose=TRUE))",
        "status <- ifelse(is(status, 'try-error'), 1, 0)",
        "quit(save='no', status=status, runLast=FALSE)",
    ]
    cmd = ["R", "--vanilla", "-e", "; ".join(R_script)]
    return subprocess.run(cmd, text=True, check=True)

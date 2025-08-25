# -*- coding: utf-8 -*-
"""Code for testing the R package installation helper function"""

from snappy_wrappers.utils import install_R_package, install_R_packages

R_packages_json = r"""[
    {
        "name": "cran",
        "repository": "cran"
    },
    {
        "name": "bioc",
        "repository": "bioconductor"
    },
    {
        "name": "package",
        "repository": "github",
        "url": "username/package/subdir@*rel.ea.se"
    },
    {
        "name": "package",
        "repository": "bitbucket",
        "url": "username/package@ref"
    },
    {
        "name": "package",
        "repository": "local",
        "url": "/path/to/package.tar.gz"
    }
]"""

def test_install_R_packages(mocker, fp, fs):
    """Tests install_R_package"""
    mocker.patch("snappy_wrappers.utils.os.makedirs", return_value=True)
    fs.create_file("/path/to/package.tar.gz")
    fs.create_file("/path/to/json", contents=R_packages_json)
    packages = [
        {
            "name": "cran",
            "repo": "cran",
            "install": "install.packages('{}', lib='/path/to/lib', repos='https://cloud.r-project.org', update=FALSE, ask=FALSE)",
            "check": "find.package('cran', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "bioc",
            "repo": "bioconductor",
            "install": "BiocManager::install('{}', lib='/path/to/lib', update=FALSE, ask=FALSE)",
            "check": "find.package('bioc', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "package",
            "url": "username/package/subdir@*rel.ea.se",
            "repo": "github",
            "install": "remotes::install_github('{}', lib='/path/to/lib', upgrade='never')",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "package",
            "url": "username/package@ref",
            "repo": "bitbucket",
            "install": "remotes::install_bitbucket('{}', lib='/path/to/lib', upgrade='never')",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "package",
            "url": "/path/to/package.tar.gz",
            "repo": "local",
            "install": "install.packages('{}', repos=NULL, lib='/path/to/lib', update=FALSE, ask=FALSE)",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
    ]
    for package in packages:
        script = "; ".join(
            [
                ".libPaths(c(.libPaths(), '/path/to/lib'))",
                package["install"].format(package.get("url", package["name"])),
                "status <- try({})".format(package["check"]),
                "status <- ifelse(is(status, 'try-error'), 1, 0)",
                "quit(save='no', status=status, runLast=FALSE)",
            ]
        )
        fp.register_subprocess(["R", "--vanilla", "-e", script], stdout="")
    install_R_packages(dest="/path/to/lib", filename="/path/to/json")

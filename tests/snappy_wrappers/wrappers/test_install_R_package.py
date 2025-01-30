# -*- coding: utf-8 -*-
"""Code for testing the R package installation helper function"""

from snappy_wrappers.utils import install_R_package


def test_install_R_package(mocker, fp):
    """Tests install_R_package"""
    mocker.patch("snappy_wrappers.utils.os.makedirs", return_value=True)
    packages = [
        {
            "name": "cran",
            "repo": "cran",
            "install": "install.packages('{}', lib='/path/to/lib', update=FALSE, ask=FALSE)",
            "check": "find.package('cran', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "bioc",
            "repo": "bioconductor",
            "install": "BiocManager::install('{}', lib='/path/to/lib', update=FALSE, ask=FALSE)",
            "check": "find.package('bioc', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "username/package/subdir@*rel.ea.se",
            "repo": "github",
            "install": "remotes::install_github('{}', lib='/path/to/lib', upgrade='never')",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "username/package@ref",
            "repo": "bitbucket",
            "install": "remotes::install_bitbucket('{}', lib='/path/to/lib', upgrade='never')",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
        {
            "name": "/path/to/package.tar.gz",
            "repo": "local",
            "install": "remotes::install_local('{}', lib='/path/to/lib', upgrade='never')",
            "check": "find.package('package', lib.loc='/path/to/lib', quiet=FALSE, verbose=TRUE)",
        },
    ]
    for package in packages:
        script = "; ".join(
            [
                ".libPaths(c(.libPaths(), '/path/to/lib'))",
                package["install"].format(package["name"]),
                "status <- try({})".format(package["check"]),
                "status <- ifelse(is(status, 'try-error'), 1, 0)",
                "quit(save='no', status=status, runLast=FALSE)",
            ]
        )
        fp.register_subprocess(["R", "--vanilla", "-e", script], stdout="")
        install_R_package(dest="/path/to/lib", name=package["name"], repo=package["repo"])

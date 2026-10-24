#!/bin/sh
# Runuran is not on conda-forge, and sequenza is installed from its repository at a fixed
# commit; install them into the environment's R library
set -eu
R --vanilla --slave << '_EOF'
lib <- R.home("library")
remotes::install_version("Runuran", version = "0.41", lib = lib, repos = "https://cloud.r-project.org", upgrade = "never")
remotes::install_bitbucket("sequenzatools/sequenza@383879d", lib = lib, upgrade = "never")
for (package in c("Runuran", "sequenza")) loadNamespace(package, lib.loc = lib)
_EOF

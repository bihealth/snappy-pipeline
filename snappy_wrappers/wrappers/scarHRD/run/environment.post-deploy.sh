#!/bin/sh
# scarHRD and the versions of copynumber and sequenza that it needs are not on CRAN,
# Bioconductor or bioconda; install them into the environment's R library
set -eu
R --vanilla --slave << '_EOF'
lib <- R.home("library")
remotes::install_github("aroneklund/copynumber@2f54c083a4a69a76fcf441a1e411085e0174b5d9", lib = lib, upgrade = "never")
remotes::install_bitbucket("sequenzatools/sequenza@3.0.0", lib = lib, upgrade = "never")
remotes::install_github("sztup/scarHRD@c98f8bc42ed0810393a98677d415114360616725", lib = lib, upgrade = "never")
for (package in c("copynumber", "sequenza", "scarHRD")) loadNamespace(package, lib.loc = lib)
_EOF

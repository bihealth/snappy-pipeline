#!/bin/bash

# Force-set TMPDIR for the BIH cluster
export OLD_TMPDIR=$TMPDIR
export TMPDIR=$HOME/scratch/tmp/$HOSTNAME/$(date +%Y%m%d)
mkdir -p $TMPDIR

# Fix the umask.
umask ug=rwx,o=

# properties = {properties}
{exec_job}

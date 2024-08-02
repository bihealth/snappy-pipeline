#!/bin/bash

# Fix the umask.
umask ug=rwx,o=

# properties = {properties}
{exec_job}

#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Read through VCF file to check for errors."""

import sys

import vcfpy

if len(sys.argv) != 2:
    print("USAGE: check_vcf.py INPUT.vcf", file=sys.stderr)
    sys.exit(1)

reader = vcfpy.Reader.from_path(sys.argv[1])
i = -1
try:
    for i, record in enumerate(reader):
        pass
except Exception as e:
    print("Error %s after %d records" % (e, i), file=sys.stderr)
    raise  # re-raise

print("Successfully read %d records from %s" % (i, sys.argv[1]), file=sys.stderr)

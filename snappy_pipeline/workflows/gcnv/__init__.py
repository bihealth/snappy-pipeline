"""This module contains the code of the gCNV integration.

The GATK gCNV caller requires a relatively complex workflow so we have extracted the
implementations from the modules implementing the workflows here.

Overall, there are the following SNAPPY related pipeline steps:

- ``helper_gcnv_build_model_wgs``
- ``helper_gcnv_build_model_targeted_seq``
- ``wgs_cnv_calling``
- ``sv_calling_targeted``

We only implement calling in CASE mode, COHORT mode is only used for building the background
model.  However, note that we run the CASE mode on all samples from a given sheet.  This may
sound paradox but the gCNV code is meant to be called on whole cohorts and there is a
considerable startup overhead in launching the CNV calling / denoising with the builtin
model.  On startup of calling in both CASE and COHORT mode, a length compilation phase is run.
"""

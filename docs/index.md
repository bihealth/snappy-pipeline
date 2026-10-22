# Snappy Pipeline Documentation

This is the documentation for the snappy pipeline.
It has five parts.

Pipeline User Docs
: Documentation for **pipeline users**.
    Start here to learn about the pipeline.
    This section starts at {ref}`quickstart`.

The somatic pipeline in detail
: Documentation for the **users of the somatic pipeline**.
    It covers adapter trimming, barcodes, the panel of normals and copy number calling.
    This section starts at {ref}`somatic_ngs`.

Pipeline Step Docs
: Documentation for the individual pipeline steps.
    Each page has a general description, the configuration, the inputs and the output files of a step.
    This section starts at {ref}`pipeline_steps_introduction`.

Pipeline Developers Docs
: Documentation for **pipeline developers**.
    Continue here if you want to fix, change or extend the pipeline.
    This section starts at {ref}`dev_intro`.

Project Info
: Housekeeping information about the project, such as the release procedure, the author list and the changelog.
    This section starts at {ref}`howto_release`.

```{note}
Even if you want to modify the pipeline, read the user documentation first.
Start at {ref}`quickstart`, or at {ref}`installation` for the developer setup.
You need to run the pipeline to test your changes.
```

```{toctree}
:caption: Pipeline Users Docs
:hidden:

quickstart
installation
usage
overview
migrations
```

```{toctree}
:caption: Somatic pipeline
:hidden:

somatic_ngs
panel_of_normals
somatic_cnv
cbioportal_export
```

```{toctree}
:caption: Pipeline Step Docs
:hidden:
:glob:

step_intro
step_generic
step/*
```

```{toctree}
:caption: Pipeline Developers Docs
:hidden:

dev_intro
dev_conventions
dev_somatic_variant_calling
dev_ngs_mapping
dev_api
```

```{toctree}
:caption: Project Info
:hidden:

contributing
howto_release
authors
history
license
```

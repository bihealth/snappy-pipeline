# Snappy Pipeline Documentation

The documentation of the snappy pipeline has five parts.

Pipeline User Docs
: For **pipeline users**.
    Read this part first.
    It starts at {ref}`quickstart`.

The somatic pipeline in detail
: For the **users of the somatic pipeline**: adapter trimming, BQSR, the panel of normals and copy number calling.
    It starts at {ref}`somatic_ngs`.

Pipeline Step Docs
: One page per pipeline step, with its description, configuration, inputs and output files.
    It starts at {ref}`pipeline_steps_introduction`.

Pipeline Developers Docs
: For **pipeline developers**, who fix, change or extend the pipeline.
    It starts at {ref}`dev_intro`.

Project Info
: About the project itself, such as the release procedure, the author list and the changelog.
    It starts at {ref}`howto_release`.

```{note}
If you want to modify the pipeline, read the user documentation first: you test your changes by running the pipeline.
Start at {ref}`quickstart`, or at {ref}`installation` for the developer setup.
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
step_overview
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

(pipeline_steps_introduction)=

# Pipeline Step Introduction

This part has one page for each pipeline step.
{ref}`generic_pipeline_step` describes what all steps have in common.

The page of a step is generated from the code of the step.
It has the description from the module docstring, the table of `depends_on` keys with the data they require, and the default configuration.
The default configuration comes from the configuration model, so it cannot differ from the code.

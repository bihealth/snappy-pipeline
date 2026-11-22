(pipeline_steps_introduction)=

# Pipeline Step Introduction

This part has one page for each pipeline step.
{ref}`generic_pipeline_step` describes what all steps have in common.

The page of a step is generated from the code of the step.
It has the step's description (its module docstring), a table of its `depends_on` keys with the data each key requires, and the default configuration.
The default configuration comes from the step's configuration model, so it always matches the code.
{ref}`step_overview` shows all steps with the data they produce and require, as a graph and a table.

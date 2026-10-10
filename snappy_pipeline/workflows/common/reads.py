"""Reads of a library as rule inputs and params (plans.md F11).

Wrappers take the FASTQ files from ``snakemake.input.reads_left`` and ``reads_right``. If another
task wrote the reads (``adapter_trimming``), only its ``out/.done`` file is a job output, so that
is the rule input and the wrappers take the file lists from ``params.args.input`` instead.
"""

from snappy_pipeline.workflows.abstract import DATA_SETS, BaseStep


def reads_are_inputs(step: BaseStep) -> bool:
    """Return whether the FASTQ files of ``step`` exist before the run, so they are inputs."""
    source = step.depends_on.reads
    return source in ("", DATA_SETS) or step.project.task(source).step == "external_data"


def _read_lists(step: BaseStep, library_name: str) -> dict[str, list[str]]:
    groups = step.read_groups(library_name)
    reads = {"reads_left": [group.left for group in groups]}
    if reads_right := [group.right for group in groups if group.right]:
        reads["reads_right"] = reads_right
    return reads


def reads_input_files(step: BaseStep, library_name: str) -> dict[str, list[str]]:
    """Return the rule inputs that provide the reads of ``library_name``."""
    if reads_are_inputs(step):
        return _read_lists(step, library_name)
    return {"reads": step.reads_input(library_name)}


def reads_params(step: BaseStep, library_name: str) -> dict[str, dict[str, list[str]]]:
    """Return the read file lists as params, or ``{}`` if the files are rule inputs."""
    if reads_are_inputs(step):
        return {}
    return {"input": _read_lists(step, library_name)}

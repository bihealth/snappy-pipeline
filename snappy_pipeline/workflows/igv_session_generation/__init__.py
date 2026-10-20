"""IGV session generation

The step writes an IGV session XML file per pedigree.
The session shows the alignments of all pedigree members and the variant calls of the index as genome tracks.
The task has no ``tool`` field.

``depends_on.alignments`` names the ``ngs_mapping`` task.
The variants come from exactly one of ``depends_on.variants`` (:ref:`step_variant_calling`), ``depends_on.annotated_variants`` (:ref:`step_variant_annotation`) or ``depends_on.phased_variants`` (:ref:`step_variant_phasing`).
The session refers to the published files in ``output/`` of these tasks by relative path.
Only children whose parents both have a DNA library get a session.
The session lists the alignments of all pedigree members with a DNA library.
Paths in the session are relative to the directory of the session file.

Outputs
-------

The task publishes ``{library_name}/out/{library_name}.igv_session.xml`` and its ``.md5`` file below ``tasks/<task>/output/``.
The path is the same as below ``work/``, and ``{library_name}`` is the DNA library of the index.
The key ``xml`` is the one downstream tasks read through ``depends_on``.
"""

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

import os
import re
from functools import cached_property

from jinja2 import Environment, FileSystemLoader
from snakemake import shell

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import IgvSessionGeneration as IgvSessionGenerationConfigModel

#: Output keys and extensions of the session file
_OUTPUTS = {"xml": ".igv_session.xml", "xml_md5": ".igv_session.xml.md5"}


class WriteIgvSessionFileStepPart(BaseStepPart):
    """Write out the IGV session XML file."""

    #: Step name
    name = "write_igv_session_file"

    #: Class available actions
    actions = ("run",)

    published_keys = tuple(_OUTPUTS)

    def __init__(self, parent):
        super().__init__(parent)
        self.previous_step = self.parent.previous_step

    def _members(self, wildcards):
        """Return the rows of the pedigree members with a DNA library."""
        index = self.parent.pedigree_index(wildcards.library_name, include_background=True)
        return self.parent.pedigree_members(index, include_background=True)

    def _get_path_bam(self, library_name):
        return self.parent.get_upstream_paths("alignments", library_name=library_name).bam

    def _get_path_vcf(self, wildcards):
        # Phased variants are stored per child, calls and annotated variants per pedigree index
        if self.previous_step == "phased_variants":
            library_name = wildcards.library_name
        else:
            library_name = self.parent.pedigree_index(
                wildcards.library_name, include_background=True
            )
        return self.parent.get_upstream_paths(self.previous_step, library_name=library_name).vcf

    @staticmethod
    def _published(path):
        """Return the published ``output/`` path of an upstream ``work/`` path.

        The session file references the published files, not the ``work/`` files other tasks read.
        """
        return re.sub(r"(^|/)work/", r"\1output/", path, count=1)

    @dictify
    def _get_input_files_run(self, wildcards):
        # BAM file of each donor
        yield "bam", [self._get_path_bam(name) for name in self._members(wildcards)["library_name"]]
        # Input file comes from previous step.
        yield "vcf", self._get_path_vcf(wildcards)

    def _get_output_files_run(self):
        return {
            key: f"work/{{library_name}}/out/{{library_name}}{ext}" for key, ext in _OUTPUTS.items()
        }

    @staticmethod
    def render_from_template(directory, template_name, **kwargs):
        loader = FileSystemLoader(directory)
        env = Environment(loader=loader)
        template = env.get_template(template_name)
        return template.render(**kwargs)

    def run(self, wildcards, output):
        """Write out the IGV session file information"""
        # Get information for writing out the template.
        members = self._members(wildcards)
        session_dir = os.path.dirname(output.xml)
        template_vars = {
            "bams": [
                {
                    "label": donor_name,
                    "path": os.path.relpath(
                        self._published(self._get_path_bam(library_name)), session_dir
                    ),
                }
                for library_name, donor_name in zip(members["library_name"], members["donor_name"])
            ],
            "vcfs": [
                {
                    "label": "Variant Calls",
                    "path": os.path.relpath(
                        self._published(self._get_path_vcf(wildcards)), session_dir
                    ),
                }
            ],
        }
        # Write output file.
        with open(output.xml, "wt") as session_file:
            xml = self.__class__.render_from_template(
                os.path.dirname(__file__), "igv_session.xml", **template_vars
            )
            print(xml, file=session_file)
        # Write MD5 sum file.
        shell(
            r"""
            pushd $(dirname {output.xml})
            md5sum $(basename {output.xml}) >$(basename {output.xml}).md5
            """.format(output=output)
        )


class IgvSessionGenerationWorkflow(BaseStep):
    """Perform IGV session generation"""

    #: Workflow name
    name = "igv_session_generation"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"igv"}))]

    config_model_class = IgvSessionGenerationConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local IGV session output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        return {"xml": f"work/{lib}/out/{lib}.igv_session.xml"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        for prev in ("phased_variants", "annotated_variants", "variants"):
            if getattr(self.config.depends_on, prev, None):
                self.previous_step = prev
                break
        else:
            raise Exception("No previous step given!")  # pragma: no cover
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((WriteIgvSessionFileStepPart,))

    @cached_property
    def output_entities(self) -> list[str]:
        """The children whose parents are both sequenced: their DNA libraries."""
        # TODO: currently only families with children are supported
        return self.children_with_both_parents()

"""Code shared by steps in ``sv_calling_{targeted,wgs}``"""

from snakemake.io import expand

from snappy_pipeline.utils import dictify, flatten, listify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.common import vcf_outputs

#: Output keys and paths of the VCF that is the result file of the SV callers, one per pedigree
VCF_OUTPUTS = vcf_outputs("work/{library_name}/out/{library_name}")


class SvCallingGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for SV calling steps"""

    #: Whether the pedigrees of background data sets get result files
    pedigrees_include_background = True

    def check_pedigree(self, pedigree: str, library_names: list[str]) -> None:
        """Check a pedigree that gets result files, such as for a warning; nothing by default"""

    @listify
    def get_result_files(self):
        """Return the outputs of the last action, the result files, for each pedigree.

        ``skip_libraries`` of the tool skips the whole pedigree of a listed library.
        """
        templates = list(flatten(self.get_output_files(self.actions[-1]).values()))
        skipped = set(self.config.get(self.name).skip_libraries)
        pedigrees = self.parent.pedigrees(self.pedigrees_include_background)
        for pedigree, library_names in pedigrees.items():
            if set(library_names) & skipped:
                continue
            self.check_pedigree(pedigree, library_names)
            for template in templates:
                yield from expand(template, library_name=pedigree)


class SvCallingGetLogFileMixin:
    """Mixin that provides the log file of SV calling steps"""

    @dictify
    def _get_log_file(self, action: str):
        """Return dict of log files in the "log" directory"""
        if hasattr(self, f"_get_log_file_infix_{action}"):
            infix = getattr(self, f"_get_log_file_infix_{action}")()
        else:
            infix = "{library_name}"
        yield "log", f"work/{infix}/log/{infix}.{action}.log"


class SvCallingWorkflow(BaseStep):
    """Base class of ``sv_calling_{targeted,wgs}``: one VCF per pedigree, named after it"""

    produces = [DataSignature(DataType.VARIANTS, frozenset({"germline", "sv"}))]
    helper_parts = (WritePedigreeStepPart,)

    #: The pedigrees of background data sets are called, too
    uses_background = True

    #: Whether the libraries are exomes and panels, whose target regions the callers take into
    #: account (Manta's workflow, MELT's ``-exome``)
    targeted: bool

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return the VCF of a pedigree and its index for downstream consumers."""
        library_name = kwargs.get("library_name", "{library_name}")
        return {
            key: VCF_OUTPUTS[key].format(library_name=library_name) for key in ("vcf", "vcf_tbi")
        }

    @classmethod
    def default_group_by(cls, config):
        """The callers write one result per pedigree, named after the pedigree."""
        return "related"

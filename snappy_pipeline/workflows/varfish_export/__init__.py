"""Export of germline variants for VarFish

The step writes the tables that ``varfish-cli`` imports into VarFish, one set per pedigree, named after the pedigree (its ``related_group``):
the pedigree file, the small variants annotated with Mehari (``mehari annotate seqvars``, then ``mehari postprocess varfish-seqvars`` for the table), with ``depends_on.structural_variants`` the structural variants (``mehari annotate strucvars``), and with ``depends_on.qc`` the BAM quality control table from the ``samtools`` reports of ``ngs_mapping`` and the target coverage of alfred.
The task has no ``tool`` field.
It reads the calls of ``variant_calling`` (``depends_on.variants``, one VCF per pedigree) and the reports of ``ngs_mapping`` (``depends_on.alignments``).
The Mehari databases come from ``depends_on.database``, an ``external_data`` task of type ``database`` tagged ``mehari``, built by mehari 0.47 or newer (mehari-data-tx 0.14 or newer); the same task can serve a ``variant_annotation`` task with ``tool: mehari``.
Its ``files`` give the paths ``hgnc_tsv`` (the HGNC table), ``transcript_db`` (one or a list), ``clinvar_db`` and ``frequency_db``.
With ``depends_on.structural_variants``, set to an :ref:`step_sv_calling_targeted` or :ref:`step_sv_calling_wgs` task, it also exports structural variants.

The ``mehari`` section configures the annotation of the small variants as for ``variant_annotation``, with the defaults of the VarFish table: one row per variant and gene with its MANE transcript (``pick: gene``, ``pick_order``), intergenic variants kept.
``mehari.assembly`` names the genome build that ``varfish-cli`` imports into.
``strucvars.rng_seed`` fixes the UUIDs of the structural variants, so that a rerun writes the same table.

Pedigrees without an index or without a DNA library of the index are skipped with a warning.
With ``depends_on.exons``, a BED file of exons, the export keeps only the small variants near the exons.

Outputs
-------

The task links these result files below ``tasks/<task>/output/``, each with a ``.md5`` file; ``{library_name}`` is the name of the pedigree:

- ``{library_name}/out/{library_name}.ped``
- ``{library_name}/out/{library_name}.seqvars.gts.tsv.gz``
- ``{library_name}/out/{library_name}.strucvars.gts.tsv.gz``, only with ``depends_on.structural_variants``
- ``{library_name}/out/{library_name}.bam_qc.tsv.gz``, only with ``depends_on.qc``, an ``ngs_data_qc`` task with ``tool: alfred``

``varfish-cli importer caseimportinfo-create`` takes these files; it ignores a database info table and the feature effects of structural variants, which the task therefore no longer writes.
"""

import warnings
from functools import cached_property

from snakemake.io import expand

from snappy_pipeline.base import SkipLibraryWarning
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import SnakemakeDict, SnakemakeDictItemsGenerator
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType
from snappy_pipeline.workflows.abstract.warnings import InconsistentPedigreeWarning
from snappy_pipeline.workflows.common.gcnv.gcnv_common import InconsistentLibraryKitsWarning

from .model import VarfishExport as VarfishExportConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: The result files of each action, with their suffixes below ``work/{library_name}/out/``
_OUTPUTS = {
    "annotate_seqvars": {"ped": ".ped", "gts": ".seqvars.gts.tsv.gz"},
    "annotate_strucvars": {"gts": ".strucvars.gts.tsv.gz"},
    "bam_qc": {"bam_qc": ".bam_qc.tsv.gz"},
}


class MehariStepPart(BaseStepPart):
    """Annotates the variants of each pedigree with Mehari and writes the VarFish tables"""

    name = "mehari"
    actions = ("annotate_seqvars", "annotate_strucvars", "bam_qc")
    wrappers = {
        "bam_qc": "mehari/bam_qc",
        "annotate_strucvars": "mehari/annotate_strucvars",
        "annotate_seqvars": "mehari/annotate_seqvars",
    }

    #: The ``strucvars`` section passes its own extra arguments
    extra_arguments_action = "annotate_seqvars"

    default_resource_usage = ResourceUsage(threads=2, runtime="1d", mem="14GB")

    result_keys = {
        action: tuple(f"{key}{md5}" for key in outputs for md5 in ("", "_md5"))
        for action, outputs in _OUTPUTS.items()
    }

    def _get_log_file(self, action: str) -> SnakemakeDict:
        return {"log": f"work/{{library_name}}/log/{{library_name}}.{action}.log"}

    def _has_result_files(self, action: str) -> bool:
        """Whether the config asks for the files of the action."""
        if action == "annotate_strucvars":
            return bool(self.parent.config.depends_on.structural_variants)
        if action == "bam_qc":
            return bool(self.parent.config.depends_on.qc)
        return True

    def get_result_files(self):
        pedigrees = self._pedigrees_with_index()
        result = []
        for action in filter(self._has_result_files, self.actions):
            for path in self.get_output_files(action).values():
                result += expand(path, library_name=pedigrees)
        return result

    def _pedigrees_with_index(self) -> list[str]:
        """Return the pedigrees whose index has a DNA library.
        The other pedigrees are skipped with a warning that names their donors."""
        result = []
        for pedigree in self.parent.pedigrees():
            if self.parent.pedigree_index(pedigree, include_background=True) is None:
                members = self.parent.pedigree_members(pedigree, include_background=True)
                names = sorted(members["donor_name"])
                msg = f"INFO: pedigree {pedigree} has no index with a DNA library (donors: {names})"
                warnings.warn(InconsistentPedigreeWarning(msg))
                continue
            result.append(pedigree)
        return result

    def _get_output_files(self, action: str) -> dict[str, str]:
        prefix = "work/{library_name}/out/{library_name}"
        files = {}
        for key, suffix in _OUTPUTS[action].items():
            files[key] = f"{prefix}{suffix}"
            files[f"{key}_md5"] = f"{prefix}{suffix}.md5"
        return files

    @staticmethod
    def _ped(wildcards) -> str:
        return "work/write_pedigree.{library_name}/out/{library_name}.ped".format(**wildcards)

    @dictify
    def _get_input_files_annotate_seqvars(self, wildcards):
        yield "ped", self._ped(wildcards)
        yield from self.parent.reference_inputs().items()

        calling = self.parent.get_upstream_paths("variants", library_name=wildcards.library_name)
        yield "vcf", calling.vcf

        database = self.parent.get_upstream_paths("database")
        yield "hgnc_tsv", database.hgnc_tsv
        yield "transcript_db", database.transcript_db
        yield "clinvar_db", database.clinvar_db
        yield "frequency_db", database.frequency_db
        if self.config.depends_on.exons:
            yield "exon_bed", self.parent.get_upstream_paths("exons").bed

    def _get_params_annotate_seqvars(self, wildcards):
        """The options of ``mehari annotate seqvars``; ``assembly`` also for the table"""
        return {"config": self.config.mehari.model_dump(by_alias=True, exclude={"extra_arguments"})}

    def _get_output_files_annotate_seqvars(self):
        return self._get_output_files("annotate_seqvars")

    @dictify
    def _get_input_files_annotate_strucvars(self, wildcards):
        yield "ped", self._ped(wildcards)
        yield from self.parent.reference_inputs().items()

        sv_config = self.parent.get_task_config("structural_variants")
        sv_caller = str(sv_config.tool)
        library_names = self.parent.pedigree_members(
            wildcards.library_name, include_background=True
        )["library_name"].tolist()

        if set(getattr(sv_config, sv_caller).skip_libraries) & set(library_names):
            msg = (
                f"Found libraries to skip in family {library_names}.  All samples will be skipped."
            )
            warnings.warn(SkipLibraryWarning(msg))
            yield "vcf", []
            return

        if sv_caller == "gcnv":
            library_kits = [
                self.parent.library_kits.get(library_name, "") for library_name in library_names
            ]
            if len(set(library_kits)) != 1:
                names_kits = list(zip(library_names, library_kits))
                msg = (
                    "Found inconsistent library kits (more than one kit!) for pedigree "
                    f"{wildcards.library_name}.  The library name/kit pairs are "
                    f"{names_kits}.  This pedigree will be SKIPPED for gcnv export (as it was "
                    "skipped for gcnv calls)."
                )
                warnings.warn(InconsistentLibraryKitsWarning(msg))
                yield "vcf", []
                return

        variants = self.parent.get_upstream_paths(
            "structural_variants", library_name=wildcards.library_name
        )
        yield "vcf", [variants.vcf]

    def _get_params_annotate_strucvars(self, wildcards):
        return {
            "assembly": self.config.mehari.assembly,
            "rng_seed": self.config.strucvars.rng_seed,
            "extra_arguments": self.config.strucvars.extra_arguments,
        }

    def _get_output_files_annotate_strucvars(self):
        return self._get_output_files("annotate_strucvars")

    @dictify
    def _get_input_files_bam_qc(self, wildcards):
        # The samtools reports of every pedigree member, and the target coverage of those with a kit
        members = self.parent.pedigree_members(wildcards.library_name, include_background=True)
        result = {"bamstats": [], "idxstats": [], "alfred_qc": []}
        for library_name in members["library_name"]:
            qc = self.parent.get_upstream_paths("alignments", library_name=library_name)
            result["bamstats"].append(qc.bamstats)
            result["idxstats"].append(qc.idxstats)
            if library_name in self.parent.coverage_qc_libraries:
                coverage = self.parent.get_upstream_paths("qc", library_name=library_name)
                result["alfred_qc"].append(coverage.alfred_qc)
        return result

    def _get_params_bam_qc(self, wildcards):
        """The library of each ``bamstats`` and ``idxstats`` input, in their order"""
        members = self.parent.pedigree_members(wildcards.library_name, include_background=True)
        return {"libraries": members["library_name"].tolist()}

    def _get_output_files_bam_qc(self) -> SnakemakeDictItemsGenerator:
        return self._get_output_files("bam_qc")


class VarfishExportWorkflow(BaseStep):
    """Perform germline variant export to VarFish"""

    name = "varfish_export"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"varfish"}))]
    #: The pedigrees include the members of background data sets
    uses_background = True
    config_model_class = VarfishExportConfigModel

    @cached_property
    def coverage_qc_libraries(self) -> set[str]:
        """The libraries that the ``depends_on.qc`` task reports on: those with a kit among its
        libraries."""
        if not (task := self.config.depends_on.qc):
            return set()
        return set(self.task_libraries(task)["library_name"]) & set(self.library_kits)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """The pedigree file and the table of the small variants of a pedigree."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/{lib}/out/{lib}"
        return {"ped": f"{prefix}.ped", "gts": f"{prefix}.seqvars.gts.tsv.gz"}

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((WritePedigreeStepPart, MehariStepPart))

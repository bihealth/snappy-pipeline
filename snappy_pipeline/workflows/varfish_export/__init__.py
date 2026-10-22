"""Export of germline variants for VarFish

The step annotates the variants of each pedigree with Mehari and writes files that ``varfish-cli`` imports into VarFish.
The task has no ``tool`` field.
It reads the quality control data of ``ngs_mapping`` (``depends_on.alignments``) and the calls of ``variant_calling`` (``depends_on.variants``).
The Mehari databases come from ``depends_on.database``, an ``external_data`` task of type ``database`` tagged ``mehari``.
Its ``files`` give the single paths ``hgnc_tsv``, ``transcript_db``, ``clinvar_db`` and ``frequency_db``.
With ``depends_on.structural_variants``, set to an :ref:`step_sv_calling_targeted` or :ref:`step_sv_calling_wgs` task, it also exports structural variants.

The step writes one result set per pedigree, named after the DNA library of the index.
Pedigrees without an index or without a DNA library of the index are skipped with a warning.

With ``depends_on.exons``, a BED file of exons, the export keeps only the small variants near the exons.

Outputs
-------

The task links these result files below ``tasks/<task>/output/``, each with a ``.md5`` file.
The paths are the same as below ``work/``.

- ``varfish_export.{library_name}/out/mehari_annotate_seqvars.{library_name}.ped``
- ``varfish_export.{library_name}/out/mehari_annotate_seqvars.{library_name}.gts.tsv.gz``
- ``varfish_export.{library_name}/out/mehari_annotate_seqvars.{library_name}.db-infos.tsv.gz``
- ``varfish_export.{library_name}/out/mehari_annotate_strucvars.{library_name}.gts.tsv.gz``
- ``varfish_export.{library_name}/out/mehari_annotate_strucvars.{library_name}.feature-effects.tsv.gz``
- ``varfish_export.{library_name}/out/mehari_annotate_strucvars.{library_name}.db-infos.tsv.gz``
- ``varfish_export.{library_name}/out/mehari_bam_qc.{library_name}.bam-qc.tsv.gz``, only if the ``ngs_mapping`` task enables ``target_coverage_report``

Here ``{library_name}`` is the DNA library of the pedigree index.
The keys ``gts`` and ``db_infos`` (small variants) are the ones downstream tasks read through ``depends_on``.
"""

import typing
import warnings

from snakemake.io import expand
from snakemake.iocontainers import Wildcards

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

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")


#: ``cohort_name`` of libraries in a pedigree whose index has no DNA library
_NO_INDEX = "unknown"


class MehariStepPart(BaseStepPart):
    """This step part is responsible for annotating the variants with Mehari"""

    name = "mehari"
    actions = ("annotate_seqvars", "annotate_strucvars", "bam_qc")
    default_resource_usage = ResourceUsage(threads=2, runtime="1d", mem="14GB")

    def _get_log_file(self, action: str) -> SnakemakeDict:
        return {
            "log": f"work/varfish_export.{{index_ngs_library}}/log/mehari_{action}.{{index_ngs_library}}.log"
        }

    def _sv_calling_is_targeted(self) -> bool:
        """Whether ``depends_on.structural_variants`` names an ``sv_calling_targeted`` task."""
        task_name = self.parent.config.depends_on.structural_variants
        return self.parent.project.task(task_name).step == "sv_calling_targeted"

    def _has_result_files(self, action: str) -> bool:
        """Whether the config asks for the files of the action."""
        if action == "annotate_strucvars":
            return bool(self.parent.config.depends_on.structural_variants)
        if action == "bam_qc":
            return self.parent.get_task_config("alignments").target_coverage_report.enabled
        return True

    def get_result_files(self):
        result = []
        for action in filter(self._has_result_files, self.actions):
            index_ngs_libraries = self._get_index_ngs_libraries(
                require_consistent_pedigree_kits=(
                    action == "annotate_strucvars" and self._sv_calling_is_targeted()
                )
            )
            for path in self.get_output_files(action).values():
                result += expand(
                    path.replace("work/", "output/", 1),
                    index_ngs_library=list(index_ngs_libraries),
                )
        return result

    @dictify
    def _get_index_ngs_libraries(
        self, *, require_consistent_pedigree_kits: bool = False
    ) -> typing.Generator[typing.Tuple[str, typing.List[str]], None, None]:
        """Return ``dict`` that maps the index DNA library name to a list of all pedigree
        member's DNA library names.
        """
        self._warn_pedigrees_without_index()
        for index in self.parent.pedigrees():
            yield (
                index,
                self.parent.pedigree_members(index, include_background=True)[
                    "library_name"
                ].tolist(),
            )

    def _warn_pedigrees_without_index(self):
        """Issue a warning for the members of pedigrees whose index has no DNA library."""
        df = self.parent.build_library_dataframe(include_background=False)
        unindexed = (df["kind"] == "germline") & (df["cohort_name"] == _NO_INDEX)
        names = sorted(df.loc[unindexed, "donor_name"])
        if names:
            msg = f"INFO: pedigree index without DNA NGS library (names: {names})"
            warnings.warn(InconsistentPedigreeWarning(msg))

    @dictify
    def _get_input_files_annotate_seqvars(self, wildcards):
        yield (
            "ped",
            "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped".format(
                **wildcards
            ),
        )
        yield from self.parent.reference_inputs().items()

        calling = self.parent.get_upstream_paths(
            "variants",
            library_name=wildcards.index_ngs_library,
        )
        yield "vcf", [calling.vcf]

        database = self.parent.get_upstream_paths("database")
        yield "hgnc_tsv", database.hgnc_tsv
        yield "transcript_db", database.transcript_db
        yield "clinvar_db", database.clinvar_db
        yield "frequency_db", database.frequency_db
        if self.config.depends_on.exons:
            yield "exon_bed", self.parent.get_upstream_paths("exons").bed

    @dictify
    def _get_output_files_annotate_seqvars(self):
        # Generate paths in "work/" directory
        prefix = (
            "work/varfish_export.{index_ngs_library}/out/"
            "mehari_annotate_seqvars.{index_ngs_library}"
        )
        work_paths = {  # annotate_seqvars will write out PED file
            "ped": f"{prefix}.ped",
            "ped_md5": f"{prefix}.ped.md5",
            "gts": f"{prefix}.gts.tsv.gz",
            "gts_md5": f"{prefix}.gts.tsv.gz.md5",
            "db_infos": f"{prefix}.db-infos.tsv.gz",
            "db_infos_md5": f"{prefix}.db-infos.tsv.gz.md5",
        }
        yield from work_paths.items()

    @dictify
    def _get_input_files_annotate_strucvars(self, wildcards):
        yield (
            "ped",
            "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped".format(
                **wildcards
            ),
        )
        yield from self.parent.reference_inputs().items()

        if not self.parent.config.depends_on.structural_variants:
            raise RuntimeError("depends_on.structural_variants is not set")
        sv_config = self.parent.get_task_config("structural_variants")
        sv_callers = [str(sv_config.tool)]
        skip_libraries = {
            sv_caller: getattr(sv_config, sv_caller).skip_libraries for sv_caller in sv_callers
        }

        library_names = self.parent.pedigree_members(
            wildcards.index_ngs_library, include_background=True
        )["library_name"].tolist()

        vcfs = []
        for sv_caller in sv_callers:
            if any(map(skip_libraries[sv_caller].__contains__, library_names)):
                msg = f"Found libraries to skip in family {library_names}.  All samples will be skipped."
                warnings.warn(SkipLibraryWarning(msg))
                continue

            if sv_caller == "gcnv":
                library_kits = [
                    self.parent.library_kits.get(library_name, "") for library_name in library_names
                ]
                if len(set(library_kits)) != 1:
                    names_kits = list(zip(library_names, library_kits))
                    msg = (
                        "Found inconsistent library kits (more than one kit!) for pedigree with "
                        f"index {wildcards.index_ngs_library}.  The library name/kit pairs are "
                        f"{names_kits}.  This pedigree will be SKIPPED for gcnv export (as it was "
                        "skipped for gcnv calls)."
                    )
                    warnings.warn(InconsistentLibraryKitsWarning(msg))
                    continue

            vcfs.append(
                self.parent.get_upstream_paths(
                    "structural_variants", library_name=wildcards.index_ngs_library
                )["vcf"]
            )
        yield "vcf", vcfs

    @dictify
    def _get_output_files_annotate_strucvars(self):
        prefix = (
            "work/varfish_export.{index_ngs_library}/out/"
            "mehari_annotate_strucvars.{index_ngs_library}"
        )
        work_paths = {
            "gts": f"{prefix}.gts.tsv.gz",
            "gts_md5": f"{prefix}.gts.tsv.gz.md5",
            "feature_effects": f"{prefix}.feature-effects.tsv.gz",
            "feature_effects_md5": f"{prefix}.feature-effects.tsv.gz.md5",
            "db_infos": f"{prefix}.db-infos.tsv.gz",
            "db_infos_md5": f"{prefix}.db-infos.tsv.gz.md5",
        }
        yield from work_paths.items()

    @dictify
    def _get_input_files_bam_qc(self, wildcards):
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "alfred_qc": []}
        for library_name in members["library_name"]:
            qc = self.parent.get_upstream_paths("alignments", library_name=library_name)
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(getattr(qc, key))
            # The coverage report covers only libraries with a kit
            if library_name in self.parent.library_kits:
                result["alfred_qc"].append(qc.alfred_qc)
        return result

    @dictify
    def _get_output_files_bam_qc(self) -> SnakemakeDictItemsGenerator:
        prefix = "work/varfish_export.{index_ngs_library}/out/mehari_bam_qc.{index_ngs_library}"
        work_paths = {
            "bam_qc": f"{prefix}.bam-qc.tsv.gz",
            "bam_qc_md5": f"{prefix}.bam-qc.tsv.gz.md5",
        }
        yield from work_paths.items()

    def _get_params_bam_qc(self, wildcards: Wildcards) -> typing.Dict[str, str]:
        """Get parameters for wrapper ``variant_annotator/bam_qc``

        Creates dictionary that links library name to identifier that should be used in output file.
        The wrapper will derive the library name from the input file name, for analysis using
        externally generated data, the values will be the sample name as provided by the external
        source (sample name). For snappy-based analysis it is redundant, both keys and values will
        be the library name.

        Dictionary expected structure:
        {
            "P001-N1-DNA1-WGS1": "P001-N1-DNA1-WGS1",
            "P002-N1-DNA1-WGS1": "P001-N1-DNA1-WGS1",
            "P003-N1-DNA1-WGS1": "P001-N1-DNA1-WGS1",
        }

        :return: Dictionary linking library name to identifier that should be used in output file.
        Key: library name; Value: identifier to be used in file.
        """
        members = self.parent.pedigree_members(wildcards.index_ngs_library, include_background=True)
        return {library_name: library_name for library_name in members["library_name"]}


class VarfishExportWorkflow(BaseStep):
    """Perform germline variant export to VarFish"""

    name = "varfish_export"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"varfish"}))]
    config_model_class = VarfishExportConfigModel

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local VarFish export output paths for downstream consumers."""
        lib = kwargs.get("library_name", "{library_name}")
        prefix = f"work/varfish_export.{lib}/out/mehari_annotate_seqvars.{lib}"
        return {
            "gts": f"{prefix}.gts.tsv.gz",
            "db_infos": f"{prefix}.db-infos.tsv.gz",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((WritePedigreeStepPart, MehariStepPart))

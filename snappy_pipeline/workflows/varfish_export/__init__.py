"""Export of germline variants for VarFish

The step annotates the variants of each pedigree with Mehari and writes files that ``varfish-cli`` imports into VarFish.
The task has no ``tool`` field.
It reads the quality control data of ``ngs_mapping`` (``depends_on.alignments``) and the calls of ``variant_calling`` (``depends_on.variants``).
With ``depends_on.structural_variants``, set to an :ref:`step_sv_calling_targeted` or :ref:`step_sv_calling_wgs` task, it also exports structural variants.

The step writes one result set per pedigree, named after the DNA library of the index.
Pedigrees without an index or without a DNA library of the index are skipped with a warning.

Outputs
-------

The task publishes these files below ``tasks/<task>/output/``, each with a ``.md5`` file.
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

import re
import typing
import warnings
from pathlib import Path

from biomedsheets.shortcuts import GermlineCaseSheet, Pedigree, is_not_background
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

# TODO: the number of restart runtimes is high because tabix in HTSJDK/Jannovar is flaky...

#: Default configuration for the somatic_variant_calling step


class MehariStepPart(BaseStepPart):
    """This step part is responsible for annotating the variants with Mehari"""

    name = "mehari"
    actions = ("annotate_seqvars", "annotate_strucvars", "bam_qc")
    default_resource_usage = ResourceUsage(threads=2, runtime="1d", mem="14GB")

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def _get_log_file(self, action: str) -> SnakemakeDict:
        return {
            "log": f"work/varfish_export.{{index_ngs_library}}/log/mehari_{action}.{{index_ngs_library}}.log"
        }

    def _sv_calling_is_targeted(self) -> bool:
        """Whether ``depends_on.structural_variants`` names an ``sv_calling_targeted`` task."""
        task_name = self.parent.config.depends_on.structural_variants
        return self.parent.project.task(task_name).step == "sv_calling_targeted"

    def _is_published(self, action: str) -> bool:
        """Whether the config asks for the files of the action."""
        if action == "annotate_strucvars":
            return bool(self.parent.config.depends_on.structural_variants)
        if action == "bam_qc":
            return self.parent.get_task_config("alignments").target_coverage_report.enabled
        return True

    def get_result_files(self):
        result = []
        for action in filter(self._is_published, self.actions):
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
        for sheet in filter(is_not_background, self.parent.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if self._is_pedigree_good(pedigree):
                    index = pedigree.index.dna_ngs_library.name
                    donors = [
                        donor.dna_ngs_library.name
                        for donor in pedigree.donors
                        if donor.dna_ngs_library
                    ]
                    yield index, donors

    def _is_pedigree_good(self, pedigree: Pedigree) -> bool:
        """Check pedigrees for inconsistencies and issue warning for any.

        :return: ``True`` if there was no inconsistency reported
        """
        msg = None
        donor_names = list(sorted(d.name for d in pedigree.donors))
        if not pedigree.index:  # pragma: no cover
            msg = f"INFO: pedigree without index (name: {donor_names})"
        elif not pedigree.index.dna_ngs_library:  # pragma: no cover
            msg = f"INFO: pedigree index without DNA NGS library (names: {donor_names})"
        if msg:
            warnings.warn(InconsistentPedigreeWarning(msg))
        return not msg

    @dictify
    def _get_input_files_annotate_seqvars(self, wildcards):
        yield (
            "ped",
            "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped".format(
                **wildcards
            ),
        )
        yield "reference", self.parent.get_upstream_paths("reference").fasta

        calling = self.parent.get_upstream_paths(
            "variants",
            library_name=wildcards.index_ngs_library,
        )
        yield "vcf", [calling.vcf]

        path_mehari_db = Path(self.config.path_mehari_db)
        prefix = path_mehari_db / self.config.release.lower()
        yield "hgnc_tsv", str(path_mehari_db / "hgnc.tsv")
        yield "transcript_db", str(prefix / "seqvars" / "txs.bin.zst")
        yield "clinvar_db", str(prefix / "seqvars" / "clinvar" / "rocksdb")
        yield "frequency_db", str(prefix / "seqvars" / "frequencies" / "rocksdb")
        if self.config.path_exon_bed:
            yield "exon_bed", self.config.path_exon_bed

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
        yield "reference", self.parent.get_upstream_paths("reference").fasta

        if not self.parent.config.depends_on.structural_variants:
            raise RuntimeError("depends_on.structural_variants is not set")
        sv_config = self.parent.get_task_config("structural_variants")
        sv_callers = [str(sv_config.tool)]
        skip_libraries = {
            sv_caller: getattr(sv_config, sv_caller).skip_libraries for sv_caller in sv_callers
        }

        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        library_names = [
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        ]

        vcfs = []
        for sv_caller in sv_callers:
            if any(map(skip_libraries[sv_caller].__contains__, library_names)):
                msg = f"Found libraries to skip in family {library_names}.  All samples will be skipped."
                warnings.warn(SkipLibraryWarning(msg))
                continue

            if sv_caller == "gcnv":
                library_kits = [
                    self.parent.ngs_library_to_kit.get(library_name, "__default__")
                    for library_name in library_names
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
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "alfred_qc": []}
        for donor in pedigree.donors:
            if not donor.dna_ngs_library:
                continue
            qc = self.parent.get_upstream_paths(
                "alignments", library_name=donor.dna_ngs_library.name
            )
            for key in result:
                result[key].append(getattr(qc, key))
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
        library_name_to_file_identifier = {}
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                library_name_to_file_identifier[donor.dna_ngs_library.name] = (
                    donor.dna_ngs_library.name
                )
        return library_name_to_file_identifier


class VarfishExportWorkflow(BaseStep):
    """Perform germline variant export to VarFish"""

    name = "varfish_export"
    produces = [DataSignature(DataType.EXPORTS, frozenset({"varfish"}))]
    config_model_class = VarfishExportConfigModel
    sheet_shortcut_class = GermlineCaseSheet

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

        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here

        # Build additional information
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        """Build mapping of NGS library to kit based on the ``ngs_mapping`` configuration"""
        cov_config = self.get_task_config("alignments").target_coverage_report
        regexes = {
            item.pattern: item.name
            for item in cov_config.path_target_interval_list_mapping
            if item.name != "__default__"
        }
        result = {}
        for sheet in self.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            if library.extra_infos.get("libraryKit"):
                                library_kit = library.extra_infos.get("libraryKit")
                                for pattern, name in regexes.items():
                                    if re.match(pattern, library_kit):
                                        yield library.name, name
        return result

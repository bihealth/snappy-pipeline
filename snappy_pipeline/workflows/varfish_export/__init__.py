"""Implementation of the ``varfish_export_step``

This step processes variant call and quality control information output of other pipeline
steps and prepares them for import into VarFish via ``varfish-cli``.

==========
Properties
==========

overall stability

    **stable**

applicable to

    germline alignment and variant calling data

generally applicable to

    short read DNA sequencing

==========
Step Input
==========

The step will read in

- quality control data from ``ngs_mapping``
- variant call output data from ``variant_calling`` and ``sv_calling_targeted``

===========
Step Output
===========

TODO

=============
Configuration
=============

By default, input from ``ngs_mapping`` and ``variant_calling`` is enabled by setting the
appropriate input paths.

You should enable ``sv_calling_targeted`` or ``sv_calling_wgs`` if you used the corresponding
steps (and have the corresponding data types).

=====================
Default Configuration
=====================

The default configuration is as follows. Note that the ``path_jannovar_ser`` parameter must
be set to point to the desired transcript annotations db as generated by ``jannovar download``.

.. include:: DEFAULT_CONFIG_varfish_export.rst
"""

import re
import typing
import warnings
from itertools import chain

from biomedsheets.shortcuts import GermlineCaseSheet, Pedigree, is_not_background
from matplotlib.cbook import flatten
from snakemake.io import Wildcards, expand

from snappy_pipeline.base import SkipLibraryWarning
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import SnakemakeDict, SnakemakeDictItemsGenerator
from snappy_pipeline.workflows.abstract.warnings import InconsistentPedigreeWarning
from snappy_pipeline.workflows.common.gcnv.gcnv_common import InconsistentLibraryKitsWarning
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_pipeline.workflows.sv_calling_targeted import SvCallingTargetedWorkflow
from snappy_pipeline.workflows.variant_calling import (
    VariantCallingGetLogFileMixin,
    VariantCallingWorkflow,
)

from .model import VarfishExport as VarfishExportConfigModel

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

# TODO: the number of restart times is high because tabix in HTSJDK/Jannovar is flaky...

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = VarfishExportConfigModel.default_config_yaml_string()


class MehariStepPart(VariantCallingGetLogFileMixin, BaseStepPart):
    """This step part is responsible for annotating the variants with Mehari"""

    name = "mehari"
    actions = ("annotate_seqvars", "annotate_strucvars", "bam_qc")

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    def get_input_files(self, action) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    def get_output_files(self, action) -> SnakemakeDict:
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def get_log_file(self, action: str) -> SnakemakeDictItemsGenerator:
        self._validate_action(action)
        prefix = (
            "work/{mapper}.varfish_export.{index_ngs_library}/log/"
            f"{{mapper}}.mehari_{action}.{{index_ngs_library}}"
        )
        key_ext = (
            ("wrapper", ".wrapper.py"),
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_params(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_params_{action}")

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
        self._validate_action(action)
        return ResourceUsage(
            threads=2,
            time="1-00:00:00",
            memory="14G",
        )

    @listify
    def get_result_files(self, action):
        # Generate templates to the output paths from action's result files.
        if action == "annotate_seqvars":
            raw_path_tpls = self._get_output_files_annotate_seqvars().values()
        elif action == "annotate_strucvars":
            # Only annotate_seqvars SVs if path to step for calling them is configured.
            if (
                not self.parent.config.path_sv_calling_targeted
                and not self.parent.config.path_sv_calling_wgs
            ):
                return
            raw_path_tpls = self._get_output_files_annotate_strucvars().values()
        elif action == "bam_qc":
            raw_path_tpls = self._get_output_files_bam_qc().values()
        # Filter the templates to the paths in the output directory.
        path_tpls = [tpl for tpl in flatten(raw_path_tpls) if tpl.startswith("output/")]

        # Create concrete paths for all pedigrees in the sample sheet.
        index_ngs_libraries = self._get_index_ngs_libraries(
            require_consistent_pedigree_kits=(
                bool(self.parent.config.path_sv_calling_targeted)
                and (action == "annotate_strucvars")
            )
        )
        kwargs = {
            "index_ngs_library": list(index_ngs_libraries.keys()),
            "mapper": [self.parent.config.tools_ngs_mapping[0]],
        }
        for path_tpl in path_tpls:
            yield from expand(path_tpl, **kwargs)

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
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

        variant_calling = self.parent.sub_workflows["variant_calling"]

        path = (
            "output/{mapper}.{var_caller}.{index_ngs_library}/out/"
            "{mapper}.{var_caller}.{index_ngs_library}.vcf.gz"
        )

        vcfs = []
        for var_caller in self.parent.config.tools_variant_calling:
            vcfs.append(
                variant_calling(path).format(
                    mapper=wildcards.mapper,
                    var_caller=var_caller,
                    index_ngs_library=wildcards.index_ngs_library,
                )
            )
        yield "vcf", vcfs

    @dictify
    def _get_output_files_annotate_seqvars(self):
        # Generate paths in "work/" directory
        prefix = (
            "work/{mapper}.varfish_export.{index_ngs_library}/out/"
            "{mapper}.mehari_annotate_seqvars.{index_ngs_library}"
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
        # Generate paths in "output/" directory
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(
                    work_paths.values(), self.get_log_file("annotate_seqvars").values()
                )
            ],
        )

    def _get_params_annotate_seqvars(self, wildcards: Wildcards) -> typing.Dict[str, typing.Any]:
        return {
            "path_exon_bed": self.config.path_exon_bed,
            "reference": self.parent.w_config.static_data_config.reference.path,
            "path_mehari_db": self.config.path_mehari_db,
        }

    @dictify
    def _get_input_files_annotate_strucvars(self, wildcards):
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"

        if self.parent.config.path_sv_calling_targeted:
            sv_calling = self.parent.sub_workflows["sv_calling_targeted"]
            sv_callers = self.parent.config.tools_sv_calling_targeted
            skip_libraries = {
                sv_caller: getattr(
                    self.parent.w_config.step_config["sv_calling_targeted"], sv_caller
                ).skip_libraries
                for sv_caller in sv_callers
            }
        elif self.parent.config.path_sv_calling_wgs:
            sv_calling = self.parent.sub_workflows["sv_calling_wgs"]
            sv_callers = self.parent.config.tools_sv_calling_wgs.dna
            skip_libraries = {
                sv_caller: getattr(
                    self.parent.w_config.step_config["sv_calling_wgs"], sv_caller
                ).skip_libraries
                for sv_caller in sv_callers
            }
        else:
            raise RuntimeError("Neither targeted nor WGS SV calling configured")

        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        library_names = [
            donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
        ]

        path = (
            "output/{mapper}.{sv_caller}.{index_ngs_library}/"
            "out/{mapper}.{sv_caller}.{index_ngs_library}.vcf.gz"
        )

        vcfs = []
        for sv_caller in sv_callers:
            if any(map(skip_libraries[sv_caller].__contains__, library_names)):
                msg = (
                    f"Found libraries to skip in family {library_names}.  All samples will be skipped "
                    f"for {sv_caller}."
                )
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
                sv_calling(path).format(
                    mapper=wildcards.mapper,
                    sv_caller=sv_caller,
                    index_ngs_library=wildcards.index_ngs_library,
                )
            )
        yield "vcf", vcfs

    @dictify
    def _get_output_files_annotate_strucvars(self):
        prefix = (
            "work/{mapper}.varfish_export.{index_ngs_library}/out/"
            "{mapper}.mehari_annotate_strucvars.{index_ngs_library}"
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
        # Generate paths in "output/" directory
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(
                    work_paths.values(), self.get_log_file("annotate_strucvars").values()
                )
            ],
        )

    #: Alias the get params function.
    _get_params_annotate_strucvars = _get_params_annotate_seqvars

    @dictify
    def _get_input_files_bam_qc(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "alfred_qc": []}
        for donor in pedigree.donors:
            if not donor.dna_ngs_library:
                continue
            tpl = (
                f"output/{wildcards.mapper}.{donor.dna_ngs_library.name}/report/bam_qc/"
                f"{wildcards.mapper}.{donor.dna_ngs_library.name}.bam.%s.txt"
            )
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(ngs_mapping(tpl % key))
            path = (
                f"output/{wildcards.mapper}.{donor.dna_ngs_library.name}/report/alfred_qc/"
                f"{wildcards.mapper}.{donor.dna_ngs_library.name}.alfred.json.gz"
            )
            result["alfred_qc"].append(ngs_mapping(path))
        return result

    @dictify
    def _get_output_files_bam_qc(self) -> SnakemakeDictItemsGenerator:
        prefix = (
            "work/{mapper}.varfish_export.{index_ngs_library}/out/"
            "{mapper}.mehari_bam_qc.{index_ngs_library}"
        )
        work_paths = {
            "bam_qc": f"{prefix}.bam-qc.tsv.gz",
            "bam_qc_md5": f"{prefix}.bam-qc.tsv.gz.md5",
        }
        yield from work_paths.items()
        yield (
            "output_links",
            [
                re.sub(r"^work/", "output/", work_path)
                for work_path in chain(work_paths.values(), self.get_log_file("bam_qc").values())
            ],
        )

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
    sheet_shortcut_class = GermlineCaseSheet

    @classmethod
    def default_config_yaml(cls):
        """Return default config YAML, to be overwritten by project-specific one"""
        return DEFAULT_CONFIG

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            config_model_class=VarfishExportConfigModel,
            previous_steps=(VariantCallingWorkflow, SvCallingTargetedWorkflow, NgsMappingWorkflow),
        )

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes((WritePedigreeStepPart, MehariStepPart, LinkOutStepPart))

        # Register sub workflows
        self.register_sub_workflow("variant_calling", self.config.path_variant_calling)
        if self.config.path_sv_calling_targeted:
            self.register_sub_workflow("sv_calling_targeted", self.config.path_sv_calling_targeted)
        if self.config.path_sv_calling_wgs:
            self.register_sub_workflow("sv_calling_wgs", self.config.path_sv_calling_wgs)
        self.register_sub_workflow("ngs_mapping", self.config.path_ngs_mapping)

        # Copy over "tools" setting from variant_calling/ngs_mapping if not set here
        step_config = self.w_config.step_config
        if not self.config.tools_ngs_mapping:
            self.config.tools_ngs_mapping = step_config["ngs_mapping"].tools.dna
        if not self.config.tools_variant_calling and "variant_calling" in step_config:
            self.config.tools_variant_calling = step_config["variant_calling"].tools
        if (
            not self.config.tools_sv_calling_targeted
            and "sv_calling_targeted" in self.w_config.step_config
        ):
            self.config.tools_sv_calling_targeted = step_config["sv_calling_targeted"].tools
        if not self.config.tools_sv_calling_wgs and "sv_calling_wgs" in self.w_config.step_config:
            self.config.tools_sv_calling_wgs = step_config["sv_calling_wgs"].tools

        # Build additional information
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()

    @dictify
    def _build_ngs_library_to_kit(self):
        """Build mapping of NGS library to kit based on the ``ngs_mapping`` configuration"""
        cov_config = self.w_config.step_config["ngs_mapping"].target_coverage_report
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

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        for action in self.sub_steps["mehari"].actions:
            yield from self.sub_steps["mehari"].get_result_files(action)

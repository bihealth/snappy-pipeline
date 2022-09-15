# -*- coding: utf-8 -*-
"""Implementation of the ``wgs_sv_export_external`` step.

The ``wgs_sv_export_external`` step takes as the input externally generated VCFs and uses
``varfish-annotator-cli annotate-sv`` commmand to create files fit for import into VarFish Server.

==========
Stability
==========

TODO

==========
Step Input
==========

TODO

===========
Step Output
===========

TODO

====================
Global Configuration
====================

TODO

=====================
Default Configuration
=====================

.. include:: DEFAULT_CONFIG_wgs_sv.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
"""

import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.base import MissingConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInExternalStepPart,
    LinkInPathGenerator,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeSampleNameStepPart,
)

#: Extension of files
EXTS = (".tsv.gz", ".tsv.gz.md5")

#: Infixes to use for file name generation
INFIXES = ("gts", "feature-effects", "db-infos")

#: VCF key to file extensions
KEY_EXT = {
    "vcf": ".vcf.gz",
    "vcf_md5": ".vcf.gz.md5",
    "tbi": ".vcf.gz.tbi",
    "tbi_md5": ".vcf.gz.tbi.md5",
}

#: Default configuration for the wgs_sv_export_external step
DEFAULT_CONFIG = r"""
# Default configuration wgs_sv_export_external.
step_config:
  wgs_sv_export_external:
    tool_ngs_mapping: null       # OPTIONAL: used to create output file prefix.
    tool_wgs_sv_calling: null    # OPTIONAL: used to create output file prefix.
    merge_vcf_flag: true         # REQUIRED: true if pedigree VCFs still need merging.
    search_paths: []             # REQUIRED: path to all VCF files.
    search_patterns: []          # REQUIRED: list of search pattern, ex.: [{"vcf": "*/*.vcf.gz"}]
    path_refseq_ser: REQUIRED    # REQUIRED: path to RefSeq .ser file
    path_ensembl_ser: REQUIRED   # REQUIRED: path to ENSEMBL .ser file
    path_db: REQUIRED            # REQUIRED: path to annotator DB file to use
"""


class VarfishAnnotatorExternalStepPart(BaseStepPart):
    """Annotate externally generated VCF file using "varfish-annotator annotate-sv"."""

    #: Step name
    name = "varfish_annotator_external"

    #: Class available actions
    actions = ("annotate", "merge_vcf")

    def __init__(self, parent):
        super().__init__(parent)
        # Build shortcut from index library name to pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)
        # Path generator for linking in
        self.path_gen = LinkInPathGenerator(
            self.parent.work_dir, self.parent.data_search_infos, self.parent.config_lookup_paths
        )
        # Define mapper+caller tag
        self.mapper_caller_tag = self._get_mapper_caller_tag()

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        self._validate_action(action)
        return getattr(self, "_get_input_files_{}".format(action))

    @listify
    def _get_input_files_merge_vcf(self, wildcards):
        """"""
        if self.config["merge_vcf_flag"]:
            pedigree = self.index_ngs_library_to_pedigree.get(wildcards.index_ngs_library)
            for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            yield f"work/input_links/{library.name}/.done"
        else:
            yield f"work/input_links/{wildcards.index_ngs_library}/.done"

    @dictify
    def _get_input_files_annotate(self, wildcards):
        """"""
        # Pedigree
        tpl = "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        yield "ped", tpl.format(**wildcards)
        # VCF
        tpl = (
            f"work/{self.mapper_caller_tag}{{index_ngs_library}}/out/"
            f"{self.mapper_caller_tag}{{index_ngs_library}}"
        )
        for key, ext in KEY_EXT.items():
            yield key, tpl.format(**wildcards) + ext

    @dictify
    def get_output_files(self, action):
        """Return output files for the export"""
        self._validate_action(action)
        return getattr(self, "_get_output_files_{}".format(action))()

    @dictify
    def _get_output_files_merge_vcf(self):
        tpl = (
            f"work/{self.mapper_caller_tag}{{index_ngs_library}}/out/"
            f"{self.mapper_caller_tag}{{index_ngs_library}}"
        )
        for key, ext in KEY_EXT.items():
            yield key, tpl + ext

    @dictify
    def _get_output_files_annotate(self):
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        for infix in INFIXES:
            key = infix.replace("-", "_")
            yield key, prefix + f".{infix}.tsv.gz"
            yield key + "_md5", prefix + f".{infix}.tsv.gz.md5"

    def get_log_file(self, action):
        self._validate_action(action)
        return getattr(self, "_get_log_file_{}".format(action))()

    @dictify
    def _get_log_file_merge_vcf(self):
        prefix = (
            f"work/{self.mapper_caller_tag}varfish_annotated.merge_vcf.{{index_ngs_library}}/log/"
            f"{self.mapper_caller_tag}varfish_annotated.merge_vcf.{{index_ngs_library}}"
        )
        key_ext = (
            ("log", ".log"),
            ("log_md5", ".log.md5"),
            ("conda_info", ".conda_info.txt"),
            ("conda_info_md5", ".conda_info.txt.md5"),
            ("conda_list", ".conda_list.txt"),
            ("conda_list_md5", ".conda_list.txt.md5"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    @dictify
    def _get_log_file_annotate(self):
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/log/"
            "varfish_annotated.{index_ngs_library}"
        )
        key_ext = (
            ("wrapper", ".wrapper.py"),
            ("wrapper_md5", ".wrapper.py.md5"),
            ("log", ".log"),
            ("log_md5", ".log.md5"),
            ("conda_info", ".conda_info.txt"),
            ("conda_info_md5", ".conda_info.txt.md5"),
            ("conda_list", ".conda_list.txt"),
            ("conda_list_md5", ".conda_list.txt.md5"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext

    def get_resource_usage(self, action):
        """Get Resource Usage

        :param action: Action (i.e., step) in the workflow, example: 'run'.
        :type action: str

        :return: Returns ResourceUsage for step.
        """
        self._validate_action(action)
        if action == "annotate":
            return ResourceUsage(
                threads=2,
                time="4-04:00:00",  # 4 days and 4 hours
                memory=f"{7 * 1024 * 2}M",
            )
        else:
            return ResourceUsage(
                threads=1,
                time="02:00:00",  # 2 hours
                memory=f"{7 * 1024 * 2}M",
            )

    def get_params(self, action):
        self._validate_action(action)
        return getattr(self, "_get_params_{}".format(action))

    def _get_params_merge_vcf(self, wildcards):
        result = {
            "input": list(sorted(self._collect_vcfs(wildcards))),
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
        }
        return result

    def _get_params_annotate(self, wildcards):
        result = {"is_wgs": True, "step_name": "wgs_sv_export_external"}
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if (
                donor.dna_ngs_library
                and donor.dna_ngs_library.extra_infos.get("libraryType") == "WGS"
            ):
                result["is_wgs"] = True
                return result
        return result

    def _collect_vcfs(self, wildcards):
        """Yield path to pedigree VCF"""
        # Seen paths list
        seen = []
        base_path_in = "work/input_links/{library_name}"
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
            folder_name = donor.dna_ngs_library.name
            for _, path_infix, filename in self.path_gen.run(
                folder_name=folder_name, pattern_set_keys=("vcf",)
            ):
                path = os.path.join(base_path_in, path_infix, filename).format(
                    library_name=donor.dna_ngs_library.name
                )
                if path in seen:
                    print(f"WARNING: ignoring path seen before {path}", file=sys.stderr)
                else:
                    seen.append(path)
                    yield path

    def _collect_sample_ids(self, wildcards):
        """Yield sample ids in pedigree"""
        # Seen paths list
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
            yield donor.name

    def _get_mapper_caller_tag(self):
        """Get mapper and caller tag

        :return: Returns tag to be used to name intermediate and final files. Tag based on
        information provided in configuration. Output examples: 'bwa.delly2.', 'dragen.', or ''.
        """
        mapper = self.config["tool_ngs_mapping"]
        caller = self.config["tool_wgs_sv_calling"]
        if mapper and caller:
            return f"{mapper}.{caller}."
        elif mapper or caller:
            mapper = mapper or ""
            caller = caller or ""
            return f"{mapper}{caller}."
        else:
            return ""


class WgsSvExportExternalWorkflow(BaseStep):
    """Perform germline WGS SV export external"""

    #: Workflow name
    name = "wgs_sv_export_external"

    #: Default biomed sheet class
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
            (),
        )
        # Load external data search information
        self.data_search_infos = list(self._load_data_search_infos())
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeSampleNameStepPart,
                VarfishAnnotatorExternalStepPart,
                LinkInExternalStepPart,
                LinkOutStepPart,
            )
        )

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        name_pattern = "varfish_annotated.{index_library.name}"
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + ".{infix}{ext}"),
            infix=INFIXES,
            ext=EXTS,
        )
        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "log", name_pattern + "{ext}"),
            ext=(
                ".log",
                ".log.md5",
                ".conda_info.txt",
                ".conda_info.txt.md5",
                ".conda_list.txt",
                ".conda_list.txt.md5",
            ),
        )

    def _yield_result_files(self, tpl, **kwargs):
        """Build output paths from path template and extension list"""
        for sheet in filter(is_not_background, self.shortcut_sheets):
            for pedigree in sheet.cohort.pedigrees:
                if not pedigree.index:  # pragma: no cover
                    msg = "INFO: pedigree without index (names: {})"
                    print(
                        msg.format(list(sorted(d.name for d in pedigree.donors))), file=sys.stderr
                    )
                    continue
                elif not pedigree.index.dna_ngs_library:  # pragma: no cover
                    msg = "INFO: pedigree index without DNA NGS library (names: {})"
                    print(
                        msg.format(  # pragma: no cover
                            list(sorted(d.name for d in pedigree.donors))
                        ),
                        file=sys.stderr,
                    )
                    continue  # pragma: no cover
                yield from expand(tpl, index_library=[pedigree.index.dna_ngs_library], **kwargs)

    def check_config(self):
        """Check configuration

        :raises: MissingConfiguration: on missing or invalid configuration.
        """
        # Initialise variables
        fail_test_bool = False
        error_msg = "Missing or invalid configuration issue(s):\n"
        required_file_keys = ("path_refseq_ser", "path_ensembl_ser", "path_db")

        # Test files
        for key in required_file_keys:
            path_ = self.config[key]
            try:
                if not os.path.isfile(path_):
                    error_msg += f"- Value for '{key}' is not a file: {path_}\n"
                    fail_test_bool = True
            except (ValueError, KeyError):
                error_msg += f"- Value '{key}' is not properly defined: {path_}\n"
                fail_test_bool = True

        # Test search paths
        search_paths = [item for item in self.config["search_paths"]]
        if len(search_paths) == 0:
            error_msg += "- Value for 'search_paths' cannot be empty.\n"
            fail_test_bool = True
        else:
            for path_ in search_paths:
                if not os.path.isdir(path_):
                    error_msg += f"- Path in 'search_paths' is not a directory: {path_}\n"
                    fail_test_bool = True

        # Test search pattern
        search_patterns = [item for item in self.config["search_patterns"]]
        if len(search_patterns) == 0:
            error_msg += "- Value for 'search_patterns' cannot be empty.\n"
            fail_test_bool = True
        else:
            for value in search_patterns:
                if not isinstance(value, dict):
                    error_msg += (
                        "- Value in 'search_patterns' is not a dictionary.\n"
                        "Expected: [{'vcf': '*/*.vcf.gz'}]\n"
                        f"Observed {type(value)}: '{value}'\n"
                    )
                    fail_test_bool = True
        # Assert
        if fail_test_bool:
            raise MissingConfiguration(error_msg)

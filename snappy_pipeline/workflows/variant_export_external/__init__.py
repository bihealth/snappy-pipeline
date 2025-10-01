# -*- coding: utf-8 -*-
"""Implementation of the ``variant_export_external`` step.

The ``variant_export_external`` takes externally generated VCF files and uses
``varfish-annotator-cli annotate`` command to create files fit for import into VarFish
Server. If externally generated BAM files are also available, it can optionally generate QC reports.

==========
Stability
==========

This step is considered experimental, use it at your own discretion.

==========
Step Input
==========

The variant export external step uses VCFs externally generated as input.

===========
Step Output
===========

For all samples, the workflow will execute ``varfish-annotator-cli``, pedigrees with multiple VCFs
will be merged based on `id`. The name of the index's primary DNA NGS library will be used as an
identification token in the output file. The following is an example of the files that
will be generated if BAM files are also available:


::
    output/
    +-- varfish_annotated.P001-N1-DNA1-WGS1
    |   `-- out
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.bam-qc.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.bam-qc.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.db-infos.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.feature-effects.tsv.gz.md5
    |   |   |-- varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz
    |   |   +-- varfish_annotated.P001-N1-DNA1-WGS1.gts.tsv.gz.md5
    |   |
    |   +-- log
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.conda_info.txt
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.conda_info.txt.md5
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.conda_list.txt
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.conda_list.txt.md5
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.log
    |       |-- varfish_annotated.annotate.P001-N1-DNA1-WGS1.log.md5
    |       |-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.conda_info.txt
    |       |-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.conda_info.txt.md5
    |       |-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.conda_list.txt
    |       |-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.conda_list.txt.md5
    |       |-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.log
    |       +-- varfish_annotated.bam_qc.P001-N1-DNA1-WGS1.log.md5
    |
    [...]

====================
Global Configuration
====================

Not applicable.

=====================
Default Configuration
=====================

.. include:: DEFAULT_CONFIG_variant_export_external.rst

==================
Parallel Execution
==================

Parallel execution is not performed currently.
"""

import os
import sys

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInBaiExternalStepPart,
    LinkInBamExternalStepPart,
    LinkInPathGenerator,
    LinkInVcfExternalStepPart,
    LinkOutStepPart,
    ResourceUsage,
    WritePedigreeSampleNameStepPart,
)
from snappy_pipeline.workflows.ngs_mapping import TargetCovReportStepPart

from .model import VariantExportExternal as VariantExportExternalConfigModel

#: Default configuration for the somatic_variant_calling step
DEFAULT_CONFIG = VariantExportExternalConfigModel.default_config_yaml_string()


class BamReportsExternalStepPart(TargetCovReportStepPart):
    """Build target coverage report and QC report for external BAM files"""

    #: Step name
    name = "bam_reports"

    #: Class available actions
    actions = ("bam_qc", "collect", "run")

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

    @staticmethod
    def _get_input_files_bam_qc(wildcards):
        yield f"work/input_links/{wildcards.library_name}/.done_bam_external"
        yield f"work/input_links/{wildcards.library_name}/.done_bai_external"

    def _get_input_files_run(self, wildcards):
        yield f"work/input_links/{wildcards.library_name}/.done_bam_external"
        yield f"work/input_links/{wildcards.library_name}/.done_bai_external"

    @dictify
    def get_output_files(self, action):
        """Return output files"""
        result = super().get_output_files(action)
        result["output_links"] = []
        yield from result.items()

    @listify
    def _get_input_files_collect(self, wildcards):
        _ = wildcards
        mapper_lib = "{mapper}.{library_name}"
        yield f"work/{mapper_lib}/report/alfred_qc/{mapper_lib}.alfred.json.gz"

    @dictify
    def _get_output_files_bam_qc_work(self):
        for report in ("bamstats", "flagstats", "idxstats"):
            report_path = (
                f"work/{{mapper}}.{{library_name}}/report/bam_qc/"
                f"{{mapper}}.{{library_name}}.bam.{report}.txt"
            )
            yield report, report_path
            yield report + "_md5", report_path + ".md5"

    def get_log_file(self, action):
        self._validate_action(action)
        if action == "run":
            return "work/{mapper}.{library_name}/log/snakemake.target_coverage.log"
        elif action == "bam_qc":
            return self._get_log_file_bam_qc()
        else:
            return "work/target_cov_report/log/snakemake.target_coverage.log"

    def get_args(self, action):
        assert action in (
            "run",
            "bam_qc",
        ), "Parameters only available for actions 'run' and 'bam_qc'."
        return getattr(self, f"_get_args_{action}")

    @staticmethod
    @dictify
    def _get_log_file_bam_qc():
        prefix = "work/{mapper}.{library_name}/log/{mapper}.{library_name}.bam_qc"
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

    def _get_args_run(self, wildcards):
        return {
            "bam": sorted(list(self._collect_bam_files(wildcards))),
            "bam_count": len(sorted(list(self._collect_bam_files(wildcards)))),
            "path_targets_bed": self.config.target_coverage_report.path_targets_bed,
        }

    def _get_args_bam_qc(self, wildcards):
        return {
            "bam": sorted(list(self._collect_bam_files(wildcards))),
            "bam_count": len(sorted(list(self._collect_bam_files(wildcards)))),
        }

    def _collect_bam_files(self, wildcards):
        """Collect BAM files."""
        # Initialise variables
        seen = []  # seen paths list
        base_path_in = "work/input_links/{library_name}"
        library_name = wildcards.library_name
        folder_name = library_name

        for _, path_infix, filename in self.path_gen.run(
            folder_name=folder_name, pattern_set_keys=("bam",)
        ):
            path = os.path.join(base_path_in, path_infix, filename).format(
                library_name=library_name
            )
            if path in seen:
                print(f"WARNING: ignoring path seen before {path}", file=sys.stderr)
            else:
                seen.append(path)
                yield path


class VarfishAnnotatorAnnotateStepPart(BaseStepPart):
    """Annotate VCF file using "varfish-annotator annotate"."""

    #: Step name
    name = "varfish_annotator_external"

    #: Class available actions
    actions = ("annotate", "bam_qc", "gvcf_to_vcf", "merge_vcf")

    #: VCF key to file extensions
    vcf_key_ext_dict = {
        "vcf": ".vcf.gz",
        "vcf_md5": ".vcf.gz.md5",
        "vcf_tbi": ".vcf.gz.tbi",
        "vcf_tbi_md5": ".vcf.gz.tbi.md5",
    }

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
        # External tool prefix - included in merged VCF name
        self.external_tool_prefix = self._get_external_source_prefix()

    def _get_external_source_prefix(self):
        """Get external source prefix

        :return: Returns external tool prefix if any provided, example: "dragen.". Otherwise,
        returns empty string.
        """
        if self.config.external_tool:
            return self.config.external_tool.lower() + "."
        return ""

    def get_input_files(self, action):
        """Return path to pedigree input file"""
        self._validate_action(action)
        return getattr(self, f"_get_input_files_{action}")

    @listify
    def _get_input_files_gvcf_to_vcf(self, wildcards):
        yield f"work/input_links/{wildcards.index_ngs_library}/.done"

    @listify
    def _get_input_files_merge_vcf(self, wildcards):
        if self.config.merge_vcf_flag:
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
        # Pedigree
        yield "ped", "work/write_pedigree.{index_ngs_library}/out/{index_ngs_library}.ped"
        # Reference
        yield "reference", self.w_config.static_data_config.reference.path
        # VCF
        tpl = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/out/"
            f"{self.external_tool_prefix}{{index_ngs_library}}"
        )
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, tpl.format(**wildcards) + ext

    def _get_input_files_bam_qc(self, wildcards):
        # Get names of primary libraries of the selected pedigree.  The pedigree is selected
        # by the primary DNA NGS library of the index.
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        result = {"bamstats": [], "flagstats": [], "idxstats": [], "cov_qc": []}
        for donor in pedigree.donors:
            mapper = self.external_tool_prefix[:-1]  # strip trailing dot
            library_name = donor.dna_ngs_library.name
            if not donor.dna_ngs_library:
                continue
            tpl = f"work/{mapper}.{library_name}/report/bam_qc/{mapper}.{library_name}.bam.%s.txt"
            for key in ("bamstats", "flagstats", "idxstats"):
                result[key].append(tpl % key)
            if donor.dna_ngs_library.name not in self.parent.ngs_library_list:
                continue
            path = f"work/{mapper}.{library_name}/report/alfred_qc/{mapper}.{library_name}.alfred.json.gz"
            result["cov_qc"].append(path)

        return result

    @dictify
    def get_output_files(self, action):
        """Return output files for action"""
        self._validate_action(action)
        return getattr(self, f"_get_output_files_{action}")()

    @dictify
    def _get_output_files_gvcf_to_vcf(self):
        tpl = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/out/"
            f"{self.external_tool_prefix}{{index_ngs_library}}"
        )
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, tpl + ext

    @dictify
    def _get_output_files_merge_vcf(self):
        tpl = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/out/"
            f"{self.external_tool_prefix}{{index_ngs_library}}"
        )
        for key, ext in self.vcf_key_ext_dict.items():
            yield key, tpl + ext

    @dictify
    def _get_output_files_annotate(self):
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        paths_work = {}
        for infix in ("gts", "db-infos"):
            key = infix.replace("-", "_")
            paths_work[key] = prefix + f".{infix}.tsv.gz"
            paths_work[key + "_md5"] = prefix + f".{infix}.tsv.gz.md5"
        paths_work["ped"] = f"{prefix}.ped"
        paths_work["ped_md5"] = f"{prefix}.ped.md5"
        yield from paths_work.items()
        yield "output_links", []

    @dictify
    def _get_output_files_bam_qc(self):
        key = "bam_qc"
        prefix = (
            "work/varfish_annotated.{index_ngs_library}/out/varfish_annotated.{index_ngs_library}"
        )
        yield key, prefix + ".bam-qc.tsv.gz"
        yield key + "_md5", prefix + ".bam-qc.tsv.gz.md5"
        yield "output_links", []

    @dictify
    def get_log_file(self, action):
        self._validate_action(action)
        if action in ("gvcf_to_vcf", "merge_vcf"):
            return self._get_log_file_complete_set(action)
        else:
            return self._get_log_file_annotation_generic(action)

    @dictify
    def _get_log_file_complete_set(self, action):
        prefix = (
            f"work/{self.external_tool_prefix}{{index_ngs_library}}/"
            f"log/{self.external_tool_prefix}{{index_ngs_library}}.{action}"
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    @staticmethod
    @dictify
    def _get_log_file_annotation_generic(action):
        prefix = (
            f"work/varfish_annotated.{{index_ngs_library}}/log/"
            f"varfish_annotated.{action}.{{index_ngs_library}}"
        )
        key_ext = (
            ("log", ".log"),
            ("conda_info", ".conda_info.txt"),
            ("conda_list", ".conda_list.txt"),
            ("wrapper", ".wrapper.py"),
            ("env_yaml", ".environment.yaml"),
        )
        for key, ext in key_ext:
            yield key, prefix + ext
            yield key + "_md5", prefix + ext + ".md5"

    def get_resource_usage(self, action: str, **kwargs) -> ResourceUsage:
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

    def get_args(self, action):
        self._validate_action(action)
        return getattr(self, f"_get_args_{action}")

    def _get_args_gvcf_to_vcf(self, wildcards):
        result = {
            "input": list(sorted(self._collect_gvcf(wildcards))),
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
            "reference_path": self.w_config.static_data_config.reference.path,
        }
        return result

    def _get_args_merge_vcf(self, wildcards):
        result = {
            "input": list(sorted(self._collect_vcfs(wildcards))),
            "sample_names": list(sorted(self._collect_sample_ids(wildcards))),
            "merge_option": self.config.merge_option,
            "reference_path": self.w_config.static_data_config.reference.path,
        }
        return result

    def _get_args_annotate(self, wildcards):
        return {"config": self.config.model_dump(by_alias=True)}

    def _get_args_bam_qc(self, wildcards):
        """Get parameters for wrapper ``variant_annotator/bam_qc``

        Creates dictionary that links library name to identifier that should be used in output file.
        The wrapper will derive the library name from the input file name, for analysis using
        externally generated data, the values will be the sample name as provided by the external
        source (sample name). For snappy-based analysis, both keys and values will be the
        library name.

        Dictionary expected structure:
        {
            "P001-N1-DNA1-WGS1": "P001",
            "P002-N1-DNA1-WGS1": "P002",
            "P003-N1-DNA1-WGS1": "P003",
        }

        :return: Dictionary linking library name to identifier that should be used in output file.
        Key: library name; Value: identifier to be used in file.
        """
        library_name_to_file_identifier = {}
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                library_name_to_file_identifier[donor.dna_ngs_library.name] = donor.secondary_id
        return library_name_to_file_identifier

    def _collect_gvcf(self, wildcards):
        """Yield path to index gVCF"""
        # Seen paths list
        seen = []
        base_path_in = "work/input_links/{library_name}"
        folder_name = wildcards.index_ngs_library
        for _, path_infix, filename in self.path_gen.run(
            folder_name=folder_name, pattern_set_keys=("vcf",)
        ):
            path = os.path.join(base_path_in, path_infix, filename).format(
                library_name=wildcards.index_ngs_library
            )
            if path in seen:
                print(f"WARNING: ignoring path seen before {path}", file=sys.stderr)
            else:
                seen.append(path)
                yield path

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
        pedigree = self.index_ngs_library_to_pedigree[wildcards.index_ngs_library]
        for donor in filter(lambda d: d.dna_ngs_library, pedigree.donors):
            yield donor.name


class VariantExportExternalWorkflow(BaseStep):
    """Perform germline variant export for externally generated VCFs"""

    #: Workflow name
    name = "variant_export_external"

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
            config_model_class=VariantExportExternalConfigModel,
            previous_steps=(),
        )
        # Load external data search information
        self.data_search_infos = list(self._load_data_search_infos())
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                WritePedigreeSampleNameStepPart,
                VarfishAnnotatorAnnotateStepPart,
                BamReportsExternalStepPart,
                LinkOutStepPart,
                LinkInVcfExternalStepPart,
                LinkInBamExternalStepPart,
                LinkInBaiExternalStepPart,
            )
        )
        self.ngs_library_list = self._build_ngs_library_list()

    @listify
    def _build_ngs_library_list(self):
        for sheet in self.shortcut_sheets:
            for donor in sheet.donors:
                for bio_sample in donor.bio_samples.values():
                    for test_sample in bio_sample.test_samples.values():
                        for library in test_sample.ngs_libraries.values():
                            yield library.name

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all primary DNA libraries and perform joint calling within pedigrees
        """
        # Initialise variables
        name_pattern = "varfish_annotated.{index.dna_ngs_library.name}"
        file_exts = (".tsv.gz", ".tsv.gz.md5")

        # Define infixes and actions - check if BAM QC is possible
        infixes = ("gts", "db-infos")
        performed_actions = ("annotate",)
        if self.config.bam_available_flag:
            infixes += ("bam-qc",)
            performed_actions += ("bam_qc",)

        yield from self._yield_result_files(
            os.path.join("output", name_pattern, "out", name_pattern + ".{infix}{ext}"),
            infix=infixes,
            ext=file_exts,
        )
        for action in performed_actions:
            name_pattern_action = name_pattern.replace(
                "varfish_annotated.", f"varfish_annotated.{action}."
            )
            prefix = os.path.join("output", name_pattern, "log", name_pattern_action + "{ext}")
            yield from self._yield_result_files(
                prefix,
                ext=(
                    ".log",
                    ".log.md5",
                    ".conda_info.txt",
                    ".conda_info.txt.md5",
                    ".conda_list.txt",
                    ".conda_list.txt.md5",
                    ".wrapper.py",
                    ".wrapper.py.md5",
                    ".environment.yaml",
                    ".environment.yaml.md5",
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
                yield from expand(tpl, index=[pedigree.index], **kwargs)

# -*- coding: utf-8 -*-
"""Export to cBioPortal

The step converts the results of the somatic tumor pipeline into the study files that cBioPortal imports.
It converts the somatic variants (``depends_on.variants``) to MAF with vcf2maf.
Optionally, it adds the copy number alterations (``copy_number_alteration``, from ``depends_on.copy_number``) and the gene expression (``expression``, from the gene counts of ``depends_on.alignments``).
The data of all samples is merged into single files per cohort.

The step has no ``tool`` field.
The task gives one result per cohort, a single study in the directory ``upload``.
Copy number data needs ``depends_on.copy_number``, and copy number and expression need ``depends_on.features``.
The segment files are written only if the study uses the reference genome ``hg19``, which follows from ``vcf2maf.ncbi_build: GRCh37``.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/upload/`` as symlinks to ``work/upload/``.

- ``meta_study.txt``, ``meta_clinical_patient.txt``, ``meta_clinical_sample.txt`` and ``meta_mutation_extended.txt``
- ``data_clinical_patient.txt``, ``data_clinical_sample.txt`` and ``data_mutation_extended.txt``
- ``case_lists/all_cases_with_mutation_data.txt``
- With ``copy_number_alteration``: ``meta_cna_gistic.txt``, ``meta_cna_log2.txt``, ``data_cna_gistic.txt``, ``data_cna_log2.txt`` and ``case_lists/all_cases_with_cna_data.txt``.
  For ``hg19``, also ``meta_segment.txt`` and ``data_segment.txt``.
  The case list ``case_lists/all_cases_with_mutation_and_cna_data.txt`` lists the samples with both.
- With ``expression``: ``meta_expression.txt``, ``data_expression.txt`` and ``case_lists/all_cases_with_mrna_rnaseq_data.txt``.
  With both options, also ``case_lists/all_cases_with_mutation_and_cna_and_mrna_data.txt``.

Downstream tasks read the keys ``meta_study`` and ``clinical_patient``.
"""

import re
import sys
from collections import namedtuple
from typing import Iterator

from biomedsheets.models import BioEntity, BioSample, NGSLibrary, SheetEntry, TestSample
from biomedsheets.shortcuts import CancerCaseSheet, CancerCaseSheetOptions, is_not_background
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, BaseStepPart, ResourceUsage
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType

from .model import CbioportalExport as CbioportalExportConfigModel

# cbioportal meta data files
META_FILES = {
    "always_present": ["meta_study.txt", "meta_clinical_patient.txt", "meta_clinical_sample.txt"],
    "sequenced": ["meta_mutation_extended.txt"],
    "cna": ["meta_cna_gistic.txt", "meta_cna_log2.txt"],
    "rna_seq_mrna": ["meta_expression.txt"],
    "segment": ["meta_segment.txt"],
}

# files holding clinical data
CLINICAL_DATA_FILES = ["data_clinical_patients.txt", "data_clinical_samples.txt"]

# files with CNA results
CNA_DATA_FILES = ["data_cna_gistic.txt", "data_cna_log2.txt", "data_segment.txt"]

# case lists files
CASE_LIST_FILES = {
    "sequenced": {
        "filename": "all_cases_with_mutation_data.txt",
        "name": "Sequenced tumors",
        "description": "Tumors with somatic variant calls",
        "stable_id": "sequenced",
        "category": "all_cases_with_mutation_data",
    },
    "cna": {
        "filename": "all_cases_with_cna_data.txt",
        "name": "Tumors with CNA data",
        "description": "Tumors with somatic Copy Number Alteration calls",
        "stable_id": "cna",
        "category": "all_cases_with_cna_data",
    },
    "rna_seq_mrna": {
        "filename": "all_cases_with_mrna_rnaseq_data.txt",
        "name": "Tumors with expression data",
        "description": "Tumors with mRNA seq expression data",
        "stable_id": "rna_seq_mrna",
        "category": "all_cases_with_mrna_rnaseq_data",
    },
    "cnaseq": {
        "filename": "all_cases_with_mutation_and_cna_data.txt",
        "name": "Sequenced tumors with CNA",
        "description": "Tumors with somatic variant & CNA calls",
        "stable_id": "cnaseq",
        "category": "all_cases_with_mutation_and_cna_data",
    },
    "3way_complete": {
        "filename": "all_cases_with_mutation_and_cna_and_mrna_data.txt",
        "name": "Sequenced tumors with CNA & expression",
        "description": "Tumors with somatic variant calls, CNA calls & expression data",
        "stable_id": "3way_complete",
        "category": "all_cases_with_mutation_and_cna_and_mrna_data",
    },
}

#: Resources of the steps that merge or convert copy number and expression tables
_TABLE_RESOURCE_USAGE = ResourceUsage(threads=2, runtime="2h", mem="8192MB")

LibraryPair = namedtuple("LibraryPair", ("normal", "tumor"))


# ================================================================================================
#
# Abstract classes: one generic, and the other for substeps dealing with a single sample at a runtime
#
# ================================================================================================


class cbioportalExportStepPart(BaseStepPart):
    """Base class for all cBioPortal export step parts

    - Builds a mapping between the tumor library id and the normal/tumor sample pair
    - Provides a list of input files based on a template, and on the iteration over the sheets.
      The template can be changed by derived classes.
    """

    #: Extraction type, must be instantiated in sub-class
    extraction_type = None

    #: input files template
    input_tpl = None

    #: Output_file name
    output_file = None

    def _yield_libraries(self):
        test_msg = 'WARNING: multiple test samples for sample "{}", test sample "{}" is ignored'
        lib_msg = 'WARNING: multiple libraries for sample "{}", library "{}" is ignored'
        samples = []
        for sheet in filter(is_not_background, self.parent.sheets):
            for donor in sheet.bio_entities.values():
                for biosample in donor.bio_samples.values():
                    if biosample.extra_infos["isTumor"]:
                        sample_name = biosample.name
                        for test_sample in biosample.test_samples.values():
                            if (
                                self.extraction_type
                                and test_sample.extra_infos["extractionType"]
                                != self.extraction_type
                            ):
                                continue
                            if sample_name in samples:
                                print(
                                    test_msg.format(sample_name, test_sample.name), file=sys.stderr
                                )
                                continue
                            samples.append(sample_name)
                            first = True
                            for lib in test_sample.ngs_libraries.values():
                                # pick only one library
                                if not first:
                                    print(lib_msg.format(sample_name, lib.name), file=sys.stderr)
                                    break
                                yield lib
                                first = False

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path of input files for merging"""
        assert self.input_tpl is not None
        for lib in self._yield_libraries():
            yield lib.test_sample.bio_sample.name, self.input_tpl.format(library_name=lib.name)

    def _samples(self):
        """Return the sample names of the merged tables, in input order"""
        return [lib.test_sample.bio_sample.name for lib in self._yield_libraries()]

    def _get_output_files_run(self):
        assert self.output_file is not None
        return self.output_file

    def _get_log_file(self, action):
        """Return path to log files for all data files"""
        return {"log": f"work/log/{self.name}.{action}.log"}

    def _get_params_run(self, wildcards):
        return dict(self.config)


class cbioportalVcf2MafStepPart(BaseStepPart):
    """Helper class for VCF2MAF step"""

    #: Step name
    name = "cbioportal_vcf2maf"

    #: Actions
    actions = ("run",)

    intermediate_outputs = {"run": ("maf",)}

    default_resource_usage = ResourceUsage(threads=2, runtime="2h", mem="5120MB")

    def __init__(self, parent):
        super().__init__(parent)

        # Build shortcut from cancer bio sample name to matched cancer sample
        donors = {}
        for lib in self.get_valid_libraries():
            ts: TestSample = lib.test_sample
            bs: BioSample = ts.bio_sample
            be: BioEntity = bs.bio_entity

            entry = donors.setdefault(be, {"tumors": [], "normal": None})

            if bs.extra_infos.get("isTumor", True):
                assert lib not in entry["tumors"], f"Duplicated library {lib.name}"
                entry["tumors"].append(lib)
            else:
                assert entry["normal"] is None, f"Multiple normals for donor {be.name}"
                entry["normal"] = lib

        self.tumor_ngs_library_to_sample_pair: dict[str, LibraryPair] = {}
        for normal_tumors in donors.values():
            normal = normal_tumors["normal"]
            for tumor in normal_tumors["tumors"]:
                self.tumor_ngs_library_to_sample_pair[tumor.name] = LibraryPair(
                    normal=normal, tumor=tumor
                )

    def get_valid_libraries(self) -> Iterator[NGSLibrary]:
        for sheet in filter(is_not_background, self.parent.sheets):
            for be in self._enabled(sheet.bio_entities):
                for bs in self._enabled(be.bio_samples):
                    for ts in self._enabled(bs.test_samples):
                        if ts.extra_infos.get("extractionType", "").upper() == "DNA":
                            for lib in self._enabled(ts.ngs_libraries):
                                yield lib

    @staticmethod
    def _enabled(items) -> Iterator[SheetEntry]:
        for i in items.values():
            if not i.disabled:
                yield i

    def _get_pair(self, wildcards: Wildcards) -> LibraryPair:
        return self.tumor_ngs_library_to_sample_pair[wildcards.tumor_library]

    def _get_output_files_run(self):
        """Return maf output file"""
        return {"maf": "work/maf/{tumor_library}/out/{tumor_library}.maf"}

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return input vcf for each output maf"""
        variants = self.parent.get_upstream_paths("variants", library_name=wildcards.tumor_library)
        yield "vcf", variants.vcf

    def _get_log_file_run(self):
        """Return path to log files for all data files"""
        return {"log": "work/maf/{tumor_library}/log/{tumor_library}.log"}

    def _get_params_run(self, wildcards):
        return {
            "tumor_sample": wildcards.tumor_library,
            "normal_sample": self._get_normal_lib_name(wildcards),
            "tumor_id": self._get_tumor_bio_sample(wildcards),
            "normal_id": self._get_normal_bio_sample(wildcards),
            "somatic_variant_annotation_tool": self.config.vcf2maf.annotation_tool,
            "ncbi_build": self.config.vcf2maf.ncbi_build,
            "Center": self.config.vcf2maf.Center,
        }

    def _get_normal_lib_name(self, wildcards):
        """Return name of normal (non-cancer) library"""
        return (n := self._get_pair(wildcards).normal) and n.name

    def _get_tumor_bio_sample(self, wildcards):
        """Return bio sample to tumor dna ngs library"""
        return self._get_pair(wildcards).tumor.test_sample.bio_sample.name

    def _get_normal_bio_sample(self, wildcards):
        """Return normal bio sample to tumor dna ngs library"""
        return (n := self._get_pair(wildcards).normal) and n.test_sample.bio_sample.name


class cbioportalMutationsStepPart(cbioportalExportStepPart):
    """Merges sample MAF files into data_mutation_extended.txt"""

    #: Step name
    name = "cbioportal_mutations"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_mutation_extended.txt"

    input_tpl = "work/maf/{library_name}/out/{library_name}.maf"


class cbioportalCns2CnaStepPart(BaseStepPart):
    """Helper class to extract gene-based log2 & copy number call from a cns file"""

    #: Step name
    name = "cbioportal_cns2cna"

    #: Actions
    actions = ("run",)

    intermediate_outputs = {"run": ("cna",)}

    default_resource_usage = _TABLE_RESOURCE_USAGE

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return the library"""
        yield "features", self.parent.get_upstream_paths("features").gtf
        copy_number = self.parent.get_upstream_paths(
            "copy_number", library_name=wildcards.tumor_library
        )
        yield "DNAcopy", copy_number.dnacopy_seg

    def _get_output_files_run(self):
        """Return cna output file"""
        return {"cna": "work/cna/{tumor_library}/out/{tumor_library}.cna"}

    def _get_log_file_run(self):
        """Return path to log files for all data files"""
        return {"log": "work/cna/{tumor_library}/log/{tumor_library}.log"}

    def _get_params_run(self, wildcards):
        return {"pipeline_id": "ENSEMBL"}


class cbioportalCnaFilesStepPart(cbioportalExportStepPart):
    """Generate cbioportal continuous (log2), discrete (pseudo-gistic) & segment tables"""

    #: Step name
    name = "cbioportal_cna"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("log2", "gistic")

    default_resource_usage = _TABLE_RESOURCE_USAGE

    input_tpl = "work/cna/{library_name}/out/{library_name}.cna"

    # Both actions merge the same per-library inputs
    @dictify
    def _get_input_files_log2(self, wildcards):
        libraries = self._yield_libraries()
        yield "tables", [self.input_tpl.format(library_name=lib.name) for lib in libraries]
        yield "mappings", self.config.path_gene_id_mappings

    _get_input_files_gistic = _get_input_files_log2

    def _get_params_log2(self, wildcards):
        return {
            "action_type": "log2",
            "samples": self._samples(),
            "extra_args": {"pipeline_id": "ENSEMBL"},
        }

    def _get_params_gistic(self, wildcards):
        return {
            "action_type": "gistic",
            "samples": self._samples(),
            "extra_args": {
                "pipeline_id": "ENSEMBL",
                "amplification": "9",
            },  # Amplification: cn >= 9 (https://doi.org/10.1038/s41586-022-04738-6)
        }

    def _get_output_files_log2(self):
        return "work/upload/data_cna_log2.txt"

    def _get_output_files_gistic(self):
        return "work/upload/data_cna_gistic.txt"


class cbioportalSegmentStepPart(cbioportalExportStepPart):
    #: Step name
    name = "cbioportal_segment"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "DNA"

    #: Actions
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_segment.txt"

    default_resource_usage = _TABLE_RESOURCE_USAGE

    def _get_input_files_run(self, wildcards):
        """Return path of input files for merging"""
        copy_number = self.parent.get_upstream_paths
        tables = [
            copy_number("copy_number", library_name=lib.name).dnacopy_seg
            for lib in self._yield_libraries()
        ]
        return {"tables": tables}

    def _get_params_run(self, wildcards):
        return {"action_type": "segment", "samples": self._samples()}


class cbioportalExpressionStepPart(cbioportalExportStepPart):
    """Generates the expression table (z-scores or FPKM)"""

    #: Step name
    name = "cbioportal_expression"

    #: Extraction type, must be instantiated in sub-class
    extraction_type = "RNA"

    #: Action
    actions = ("run",)

    #: Output file name
    output_file = "work/upload/data_expression.txt"

    default_resource_usage = _TABLE_RESOURCE_USAGE

    @dictify
    def _get_input_files_run(self, wildcards):
        """Return path of input files for merging"""
        alignments = self.parent.get_upstream_paths
        yield (
            "tables",
            [
                alignments("alignments", library_name=lib.name).gene_counts
                for lib in self._yield_libraries()
            ],
        )
        yield "mappings", self.config.path_gene_id_mappings
        yield "features", self.parent.get_upstream_paths("features").gtf

    def _get_params_run(self, wildcards):
        return {
            "action_type": "expression",
            "samples": self._samples(),
            "extra_args": {"pipeline_id": "ENSEMBL"},
        }


class cbioportalMetaFilesStepPart(BaseStepPart):
    """Generate cbioportal meta data files"""

    #: Step name
    name = "cbioportal_meta_files"

    #: Actions
    actions = ("run",)

    _get_log_file = cbioportalExportStepPart._get_log_file

    @listify
    def _get_output_files_run(self):
        files = [*META_FILES["always_present"], *META_FILES["sequenced"]]
        if self.config.copy_number_alteration.enabled:
            files += META_FILES["cna"]
            if self.config.study.reference_genome == "hg19":
                files += META_FILES["segment"]
        if self.config.expression.enabled:
            files += META_FILES["rna_seq_mrna"]
        for f in files:
            yield f"work/upload/{f}"

    def _get_params_run(self, wildcards):
        return self.config.study.model_dump(by_alias=True)


class cbioportalClinicalDataStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_clinical_data"

    #: Actions
    actions = ("run",)

    def _get_params_run(self, wildcards):
        donors = {}
        for extraction_type in ("DNA", "RNA"):
            if extraction_type == "RNA" and not self.config.expression.enabled:
                continue
            self.extraction_type = extraction_type
            for lib in self._yield_libraries():
                assert lib.test_sample.extra_infos["extractionType"] == extraction_type
                donor_name = lib.test_sample.bio_sample.bio_entity.name
                sample_name = lib.test_sample.bio_sample.name
                if donor_name not in donors:
                    donors[donor_name] = {}
                if sample_name not in donors[donor_name]:
                    donors[donor_name][sample_name] = {}
                # Multiple libraries should not be returned by _yield_libraries
                assert extraction_type not in donors[donor_name][sample_name]
                donors[donor_name][sample_name][extraction_type] = lib.name
        config = self.config.model_dump(by_alias=True, exclude={"path_gene_id_mappings"})
        return {"donors": donors, "config": config}

    def _get_output_files_run(self):
        return {
            data_type: f"work/upload/data_clinical_{data_type}.txt"
            for data_type in ("patient", "sample")
        }


class cbioportalCaseListsStepPart(cbioportalExportStepPart):
    """Generate cbioportal patient data file"""

    #: Step name
    name = "cbioportal_case_lists"

    #: Actions
    actions = ("run",)

    def _get_params_run(self, wildcards):
        samples = dict(zip(CASE_LIST_FILES.keys(), [[] for _ in range(len(CASE_LIST_FILES))]))
        self.extraction_type = "DNA"
        for lib in self._yield_libraries():
            sample_name = lib.test_sample.bio_sample.name
            samples["sequenced"] += [sample_name]
            if self.config.copy_number_alteration.enabled:
                samples["cna"] += [sample_name]
                samples["cnaseq"] += [sample_name]
        if self.config.expression.enabled:
            self.extraction_type = "RNA"
            for lib in self._yield_libraries():
                sample_name = lib.test_sample.bio_sample.name
                samples["rna_seq_mrna"] += [sample_name]
        args = {}
        for k, v in samples.items():
            if len(v) > 0:
                args[k] = CASE_LIST_FILES[k]
                args[k]["samples"] = v
        if self.config.copy_number_alteration.enabled and self.config.expression.enabled:
            args["3way_complete"] = CASE_LIST_FILES["3way_complete"]
            args["3way_complete"]["samples"] = []
            for sample_name in args["cnaseq"]["samples"]:
                if sample_name in args["rna_seq_mrna"]["samples"]:
                    args["3way_complete"]["samples"] += [sample_name]
        args["__cancer_study_id"] = self.config.study.cancer_study_id
        return args

    @dictify
    def _get_output_files_run(self):
        case_lists = {}
        case_lists["sequenced"] = CASE_LIST_FILES["sequenced"]["filename"]
        if self.config.copy_number_alteration.enabled:
            case_lists["cna"] = CASE_LIST_FILES["cna"]["filename"]
            case_lists["cnaseq"] = CASE_LIST_FILES["cnaseq"]["filename"]
        if self.config.expression.enabled:
            case_lists["rna_seq_mrna"] = CASE_LIST_FILES["rna_seq_mrna"]["filename"]
            if self.config.copy_number_alteration.enabled:
                case_lists["3way_complete"] = CASE_LIST_FILES["3way_complete"]["filename"]
        for case, filename in case_lists.items():
            yield case, f"work/upload/case_lists/{filename}"


# ================================================================================================
#
# cBioPortal step class
#
# ================================================================================================


class cbioportalExportWorkflow(BaseStep):
    """Perform cbioportal preparation"""

    #: Workflow name
    name = "cbioportal_export"

    config_model_class = CbioportalExportConfigModel

    produces = [DataSignature(DataType.EXPORTS, frozenset({"cbioportal"}))]

    #: Default biomed sheet class
    sheet_shortcut_class = CancerCaseSheet

    sheet_shortcut_kwargs = {
        "options": CancerCaseSheetOptions(allow_missing_normal=True, allow_missing_tumor=True)
    }

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local cBioPortal export output paths for downstream consumers."""
        _ = kwargs
        return {
            "meta_study": "work/upload/meta_study.txt",
            "clinical_patient": "work/upload/data_clinical_patient.txt",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)

        # cBioPortal requires the genome release as GRC[hm]3[78] in the MAF file
        # and hg19 or hg38 in the meta study & meta segment files.
        # Note that segmentation visualisation is not yet implemented for hg38 (2023-06-21)
        # For the user's convenience, the configuration is augmented automatically,
        # before the sub steps are registered, so they all have the updated config.
        translated = "unknown"
        if self.config.vcf2maf.ncbi_build == "GRCh37":
            translated = "hg19"
        if self.config.vcf2maf.ncbi_build == "GRCh38":
            translated = "hg38"
        if self.config.vcf2maf.ncbi_build in ("mm9", "mm10", "GRCm37", "GRCm38", "GRCm39"):
            translated = "mouse"
        self.config.study.reference_genome = translated

        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (
                cbioportalMetaFilesStepPart,
                cbioportalCaseListsStepPart,
                cbioportalClinicalDataStepPart,
                cbioportalCns2CnaStepPart,
                cbioportalCnaFilesStepPart,
                cbioportalSegmentStepPart,
                cbioportalMutationsStepPart,
                cbioportalVcf2MafStepPart,
                cbioportalExpressionStepPart,
            )
        )
        # Sub-workflows were registered above.

    @listify
    def get_result_files(self):
        """Return the upload files, published below ``output/``"""
        for path in self._get_work_files():
            yield re.sub(r"^work/", "output/", path)

    def _get_work_files(self):
        yield from self.sub_steps["cbioportal_meta_files"].get_output_files("run")
        yield from self.sub_steps["cbioportal_clinical_data"].get_output_files("run").values()
        yield from self.sub_steps["cbioportal_case_lists"].get_output_files("run").values()

        result_files = []
        result_files += [self.sub_steps["cbioportal_mutations"].get_output_files("run")]
        if self.config.copy_number_alteration.enabled:
            result_files += [self.sub_steps["cbioportal_cna"].get_output_files("log2")]
            result_files += [self.sub_steps["cbioportal_cna"].get_output_files("gistic")]
            if self.config.study.reference_genome == "hg19":
                result_files += [self.sub_steps["cbioportal_segment"].get_output_files("run")]
        if self.config.expression.enabled:
            result_files += [self.sub_steps["cbioportal_expression"].get_output_files("run")]

        yield from result_files

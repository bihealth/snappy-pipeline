"""Read mapping of DNA and RNA libraries

The step aligns the reads of each library to the reference genome and publishes sorted, indexed BAM files.
It also computes alignment statistics and, if configured, target coverage and fingerprints.
Marking of duplicates and trimming of adapters depend on the tool configuration.

The ``tool`` field selects the aligner:

- ``bwa`` and ``bwa_mem2``: short-read DNA.
- ``minimap2``: long-read DNA (PacBio, Oxford Nanopore).
- ``star``: short-read RNA. It also publishes gene counts and the strandedness decision.
- ``mbcs``: short-read DNA with molecular barcodes and optional base quality score recalibration (BQSR). See the notes below.

The task maps all libraries of its sample sheets, including those of background data sets, unless ``library_selection`` narrows them.
A task with a DNA tool maps the DNA libraries, a task with ``star`` the RNA libraries.

Reads come from ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
The mapper index comes from ``depends_on.index``, for example an ``external_data`` task or a ``reference_index`` task (:ref:`step_reference_index`).
A ``reference_index`` task must build it from the genome in ``depends_on.reference``.

Outputs
-------

The task publishes the following files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.bam`` and ``{library_name}/out/{library_name}.bam.bai``
- ``{library_name}/report/bam_qc/{library_name}.bam.bamstats.txt``, ``.bam.flagstats.txt`` and ``.bam.idxstats.txt``
- ``{library_name}/report/alfred_qc/{library_name}.alfred.json.gz``, if ``target_coverage_report`` is enabled and the library has a kit
- ``{library_name}/report/cov/{library_name}.cov.vcf.gz`` (with ``.tbi``), ``.cov.bw`` and ``.mq.bw``, if ``bam_collect_doc`` is enabled
- ``{library_name}/report/fingerprint/{library_name}.npz``, if ``ngs_chew_fingerprint`` is enabled
- with ``star``, for RNA libraries: ``{library_name}/out/{library_name}.GeneCounts.tab``, ``{library_name}/out/{library_name}.Junctions.tab`` and ``{library_name}/strandedness/{library_name}.decision.json``

Downstream tasks read the keys ``bam``, ``bai``, ``bamstats``, ``flagstats``, ``idxstats`` and ``alfred_qc``, and for RNA libraries ``gene_counts`` and ``strandedness``.

Reads and sample sheet
----------------------

A data set or an ``external_data`` task locates the reads below its ``search_paths``, in a folder per library.
Each entry of ``search_patterns`` holds regular expressions for the path below that folder.
The keys ``left`` and ``right`` select the mates, and a ``(?P<readgroup>...)`` group pairs the mates of one read group.

By default, every pattern must match as many ``left`` as ``right`` files.
With ``mixed_se_pe: true`` in the data set or the ``external_data`` task, ``right`` may match nothing.
Then a library is either single-end or paired-end, and different libraries may differ.

Exome enrichment kits
~~~~~~~~~~~~~~~~~~~~~

The target coverage report and the AGeNT barcode processing of ``mbcs`` read the region files of the library's enrichment kit.
The files come from an ``external_data`` task in ``kits`` mode, named in ``depends_on.kits``.
The sample sheet gets the optional library column ``libraryKit`` with the name of the kit in that task:

::

    [Custom Fields]
    key     annotatedEntity docs    type    minimum maximum unit    choices pattern
    libraryKit      ngsLibrary      exome enrichment kit    string  0       0       0       0       0

    [Data]
    patientName     sampleName      isTumor extractionType  libraryType     folderName      libraryKit
    case001 N1      N       DNA     WES     case001-N1-DNA1-WES1    V5_noUTR
    case002 N1      N       DNA     WES     case002-N1-DNA1-WES1    V5_UTR
    case003 N1      N       DNA     WES     case003-N1-DNA1-WES1
    [...]

A library uses the kit named exactly like its ``libraryKit``.
The target coverage report reads the ``targets`` file of the kit and needs ``depends_on.kits``.
Only targeted libraries (library type ``WES`` or ``Panel_seq``) have a kit.
The report skips all others.
A ``libraryKit`` that the kits task lacks is an error.
``mbcs`` with ``use_barcodes`` reads the ``baits`` file, and every DNA library needs a kit, so it fails for libraries that are not targeted.

.. code-block:: yaml

    - step: external_data
      name: kits
      config:
        produces: {type: regions, tags: [kits]}
        kits:
          V5_noUTR:
            targets: resources/v5.bed
            baits: resources/v5_baits.bed
          V5_UTR:
            targets: resources/v5_utr.bed
            baits: resources/v5_utr_baits.bed
    - step: ngs_mapping
      name: mapping
      config:
        depends_on:
          kits: kits
        target_coverage_report:
          enabled: true

Notes
-----

STAR
~~~~

STAR needs a gene model for the gene counts and for the transcriptome alignment. It comes from the ``features`` task in ``depends_on``.
STAR writes counts for unstranded, forward and reverse protocols.
The step keeps the column of the library's strand and publishes it as ``GeneCounts.tab``.
The strand is the configured protocol (0 unstranded, 1 forward, 2 reverse).
Without it, the step infers the strand with ``infer_experiment`` from RSeQC.
``infer_experiment`` runs in both cases and reads the BED file of ``depends_on.exons``: non-overlapping, sorted intervals of ubiquitously expressed genes with well-defined strand (12 columns, 6 also work).

With ``transcriptome: true``, STAR also writes the reads mapped to the transcriptome to a BAM file.
It is not indexed, because it has no positional order.

The ``mbcs`` tool
~~~~~~~~~~~~~~~~~

``mbcs`` maps with one of the DNA mappers named in ``mbcs.mapping_tool`` (``bwa`` or ``bwa_mem2``).
Handling of molecular barcodes (``use_barcodes``) and BQSR (``recalibrate``) is optional.
The mapper section (``bwa`` or ``bwa_mem2``) configures the mapping.
Use ``extra_arguments: ["-C"]`` there when barcodes are present and processed by AGeNT.

Only the proprietary `AGeNT <https://www.agilent.com/en/product/next-generation-sequencing/hybridization-based-next-generation-sequencing-ngs/ngs-software/agent-232879>`_ is implemented for barcodes.
It is not open source and not in ``bioconda``, so it must be installed separately.
The AGeNT baits file is the ``baits`` file of the library's kit in ``depends_on.kits``.
BQSR (``recalibrate``, the default) needs ``depends_on.common_variants``: a VCF of common germline variants, for example ``small_exac_common_3.vcf`` of the GATK resource bundle.

:ref:`somatic_ngs` has a complete example.

A Snakefile runs the sub-steps (barcode trimming, mapping, consensus merging, BQSR) inside one job.
The other tools use shell script wrappers.
This avoids long shell scripts and a shared temporary directory between jobs.

Fingerprints
~~~~~~~~~~~~

``ngs_chew_fingerprint`` writes a fingerprint of carefully selected germline variants with `ngs-chew <https://github.com/bihealth/ngs-chew>`_.
Fingerprints detect sample swaps and contamination.
They assign a sample to a patient more sensitively than HLA typing, which is matched between donor and recipient in stem-cell treatments.
"""

import os
import re
from functools import cached_property
from typing import Any

from snakemake.io import expand
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import InvalidConfiguration
from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    DataSignature,
    ResourceUsage,
    select_libraries,
)
from snappy_pipeline.workflows.abstract.protocol import DataType
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

from .model import NgsMapping as NgsMappingConfigModel
from .model import Tool

# TODO: Need something smarter still for @RG

#: Available read mappers for (short/Illumina) DNA-seq data
READ_MAPPERS_DNA = ("bwa", "bwa_mem2", "mbcs")

#: Available read mappers for (short/Illumina) RNA-seq data
READ_MAPPERS_RNA = ("star",)

#: Available read mappers for (long/PacBio/Nanopoare) DNA-seq data
READ_MAPPERS_DNA_LONG = ("minimap2",)

#: Files of a BWA index below its prefix. The wrapper takes the prefix from the first.
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")

#: Files of a BWA-MEM2 index below its prefix. The wrapper takes the prefix from the first.
BWA_MEM2_INDEX_EXTENSIONS = (".amb", ".ann", ".0123", ".bwt.2bit.64", ".pac")

#: Files of a STAR index directory that mark it as the index
STAR_INDEX_FILES = ("Genome", "SA", "SAindex")


class MappingGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for mapping steps"""

    tool_category = None

    #: Keys of the output files that are not published
    unpublished_keys: tuple[str, ...] = ()

    def skip_result_files_for_library(self, library_name: str) -> bool:
        """Override function such that the mapper is applicable to the library"""
        library_tool_category = self.parent.library_tool_category(library_name)
        if self.tool_category not in ("__any__", library_tool_category):
            return True
        else:
            return self.name != self.config.tool

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return the paths returned by ``self.get_output_files()`` for all
        actions in ``self.actions``, with the prefix ``work/`` replaced by ``output/``. Keys in
        ``unpublished_keys`` are skipped.

        You can override the ``skip_result_files_for_library()`` method to skip result files for
        a library.
        """
        # Skip if step part has a tool category and it is not enabled
        if self.tool_category != "__any__" and self.name != self.config.tool:
            return

        for action in self.actions:
            # Get list of result path templates.
            result_paths_tpls = [
                re.sub(r"^work/", "output/", path)
                for key, path in self.get_output_files(action).items()
                if key not in self.unpublished_keys
            ]
            #: Generate all concrete output paths.
            for path_tpl in result_paths_tpls:
                for library_name in self.parent.libraries["library_name"]:
                    if not self.skip_result_files_for_library(library_name):
                        yield from expand(path_tpl, library_name=library_name)


class ReportGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for report generation steps"""

    def skip_result_files_for_library(self, library_name: str) -> bool:
        if not getattr(self, "tool_categories", None):
            return False

        library_tool_category = self.parent.library_tool_category(library_name)
        return library_tool_category not in self.tool_categories

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return the paths returned by ``self.get_output_files()`` for all
        actions in ``self.actions``, with the prefix ``work/`` replaced by ``output/``.

        You can override the ``skip_result_files_for_library()`` method to skip result files for
        a library.
        """
        # Skip if step part has a tool category and it is not enabled
        for action in self.actions:
            # Get list of result path templates.
            result_paths_tpls = [
                re.sub(r"^work/", "output/", path)
                for path in self.get_output_files(action).values()
            ]
            #: Generate all concrete output paths.
            for library_name in self.parent.libraries["library_name"]:
                for sub_step in self.parent.sub_steps.values():
                    if (
                        isinstance(sub_step, ReadMappingStepPart)
                        and not sub_step.skip_result_files_for_library(library_name)
                        and not self.skip_result_files_for_library(library_name)
                    ):
                        for path_tpl in result_paths_tpls:
                            yield path_tpl.format(library_name=library_name)
                if action == "collect":
                    break  # only once


class ReadMappingStepPart(MappingGetResultFilesMixin, BaseStepPart):
    """Base class for read mapping step parts"""

    #: Class available actions. Read mapping step parts only support action "run".
    actions = ("run",)

    base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "sample_name": wildcards.library_name,
            "platform": "ILLUMINA",
        } | reads_params(self.parent, wildcards.library_name)

    def _get_input_files_run(self, wildcards):
        return reads_input_files(self.parent, wildcards.library_name)

    @dictify
    def _get_output_files_run(self):
        """Return the output files that all read mapping sub steps write."""
        if self.name != self.config.tool:
            return
        for ext in (".bam", ".bam.bai", ".bam.md5", ".bam.bai.md5"):
            yield ext[1:].replace(".", "_"), self.base_path_out.format(ext=ext)
        reports = ("bamstats", "flagstats", "idxstats")
        for suffix in ("", ".md5"):
            for report in reports:
                key = f"report_{report}_txt" + suffix.replace(".", "_")
                yield (
                    key,
                    f"work/{{library_name}}/report/bam_qc/{{library_name}}.bam.{report}.txt{suffix}",
                )

    def _get_log_file_run(self):
        if self.name != self.config.tool:
            return {}
        return {"log": "work/{library_name}/log/{library_name}.mapping.log"}


class BwaStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using BWA"""

    #: Step name
    name = "bwa"

    #: Tool category
    tool_category = "dna"

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.bwa.num_threads_align
        return ResourceUsage(threads=threads, runtime="3d", mem=f"{int(4.5 * 1024 * threads)}MB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        parent_args = super()._get_params_run(wildcards)
        parent_args.update(self.config.bwa.model_dump(by_alias=True))
        return parent_args

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {"index": self.parent.get_index_files("bwa")}


class BwaMem2StepPart(ReadMappingStepPart):
    name = "bwa_mem2"
    tool_category = "dna"

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.bwa_mem2.num_threads_align
        return ResourceUsage(threads=threads, runtime="3d", mem=f"{int(4.5 * 1024 * threads)}MB")

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        parent_args = super()._get_params_run(wildcards)
        parent_args.update(self.config.bwa_mem2.model_dump(by_alias=True))
        return parent_args

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {"index": self.parent.get_index_files("bwa_mem2")}


class MBCsStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment on MBC data"""

    name = "mbcs"
    tool_category = "dna"

    default_resource_usage = ResourceUsage(threads=1, runtime="72h", mem="4GB", partition="medium")

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards) | {
            **self.parent.reference_inputs(),
            "index": self.parent.get_index_files(str(self.config.mbcs.mapping_tool)),
        }
        if self.config.mbcs.use_barcodes:
            agent = getattr(self.config, self.config.mbcs.barcode_tool)
            inputs["agent_prepare"] = agent.prepare.path
            inputs["agent_mark_duplicates"] = agent.mark_duplicates.path
            inputs["agent_baits"] = self.parent.library_kit_paths(
                "kits", wildcards.library_name
            ).baits
        if self.config.mbcs.recalibrate:
            inputs |= self.parent.vcf_inputs("common_variants")
        return inputs

    def _get_params_run(self, wildcards: Wildcards):
        args = super()._get_params_run(wildcards)
        mapper = getattr(self.config, self.config.mbcs.mapping_tool)
        args |= {
            "config": self.config.mbcs.model_dump(by_alias=True),
            "mapper_config": mapper.model_dump(by_alias=True),
        }
        if self.config.mbcs.use_barcodes:
            args["barcode_config"] = getattr(self.config, self.config.mbcs.barcode_tool).model_dump(
                by_alias=True,
                exclude={"prepare": {"path"}, "mark_duplicates": {"path"}},
            )
        return args


class StarStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using STAR"""

    #: Step name
    name = "star"

    #: Tool category
    tool_category = "rna"

    #: The gene counts of STAR (all strands) are input of the strandedness step, which publishes
    #: the counts of the inferred strand
    unpublished_keys = ("gene_counts",)

    intermediate_outputs = {"run": ("gene_counts",)}

    def _get_output_files_run(self):
        """Override base class' function to make Snakemake aware of extra files for STAR."""
        output_files = super()._get_output_files_run()
        if (cfg := self.config.get(self.name)) is None:
            return output_files

        output_files["gene_counts"] = self.base_path_out.format(ext=".ReadsPerGene.tab")
        output_files["junctions"] = self.base_path_out.format(ext=".Junctions.tab")
        output_files["junctions_md5"] = output_files["junctions"] + ".md5"
        if cfg.transcriptome:
            output_files["transcriptome"] = self.base_path_out.format(ext=".toTranscriptome.bam")
            output_files["transcriptome_md5"] = output_files["transcriptome"] + ".md5"
        return output_files

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        parent_args = super()._get_params_run(wildcards)
        parent_args.update(self.config.star.model_dump(by_alias=True))
        return parent_args

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {
            "features": self.parent.get_upstream_paths("features").gtf,
            "index": self.parent.get_index_files("star"),
        }

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.star.num_threads_align
        return ResourceUsage(threads=threads, runtime="2d", mem=f"{int(3.5 * threads)}GB")


class StrandednessStepPart(BaseStepPart):
    name = "strandedness"
    actions = ("infer", "counts")

    def _get_input_files_infer(self, wildcards):
        # Must reference the mapper-specific work directory
        inputs = {"bam": "work/{library_name}/out/{library_name}.bam".format(**wildcards)}
        inputs["exon_bed"] = self.parent.get_upstream_paths("exons").bed
        return inputs

    def _get_input_files_counts(self, wildcards):
        return {
            "counts": "work/{library_name}/out/{library_name}.ReadsPerGene.tab".format(**wildcards),
            "decision": "work/{library_name}/strandedness/{library_name}.decision.json".format(
                **wildcards
            ),
        }

    def _get_output_files_infer(self):
        tpl = "work/{library_name}/strandedness/{library_name}"
        return {
            "tsv": tpl + ".infer.txt",
            "decision": tpl + ".decision.json",
            "decision_md5": tpl + ".decision.json.md5",
        }

    def _get_output_files_counts(self):
        return {
            "counts": "work/{library_name}/out/{library_name}.GeneCounts.tab",
            "counts_md5": "work/{library_name}/out/{library_name}.GeneCounts.tab.md5",
        }

    def get_result_files(self):
        if self.config.tool != "star":
            return
        tpl_out = "output/{library_name}/out/{library_name}.GeneCounts.tab"
        tpl_strandedness = "output/{library_name}/strandedness/{library_name}.decision.json"
        for library_name in self.parent.libraries["library_name"]:
            info = self.parent.library_infos[library_name]
            if info["extraction_type"] == "rna":
                yield tpl_out.format(library_name=library_name)
                yield tpl_out.format(library_name=library_name) + ".md5"
                yield tpl_strandedness.format(library_name=library_name)
                yield tpl_strandedness.format(library_name=library_name) + ".md5"

    def _get_log_file(self, action):
        return {"log": "work/{library_name}/log/{library_name}.strandedness.log"}

    def _get_params_infer(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "config": self.config.strandedness.model_dump(by_alias=True),
            "library_name": wildcards.library_name,
        }


class Minimap2StepPart(ReadMappingStepPart):
    """Support for performing long-read alignment using minimap2"""

    #: Step name
    name = "minimap2"

    #: Tool category
    tool_category = "dna_long"

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.minimap2.mapping_threads
        return ResourceUsage(threads=threads, runtime="2d", mem=f"{int(3.5 * threads)}GB")

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {"index": self.parent.get_index_files("minimap2")}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        params = super()._get_params_run(wildcards)
        params |= self.config.minimap2.model_dump(by_alias=True)
        params["library_kit"] = self.parent.library_infos[wildcards.library_name]["library_kit"]
        params["library_name"] = wildcards.library_name
        return params


class TargetCovReportStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Build target coverage report"""

    #: Step name
    name = "target_coverage_report"

    #: Run for "dna" and "dna_long" only
    tool_categories = ("dna", "dna_long")

    #: Class available actions
    actions = ("run",)

    def __init__(self, parent):
        super().__init__(parent)

    def skip_result_files_for_library(self, library_name: str) -> bool:
        return (
            not self.config.target_coverage_report
            or not self.config.target_coverage_report.enabled
            or library_name not in self.parent.library_kits
            or super().skip_result_files_for_library(library_name)
        )

    @dictify
    def _get_input_files_run(self, wildcards):
        library_name = wildcards.library_name
        yield "bam", f"work/{library_name}/out/{library_name}.bam"
        yield "bai", f"work/{library_name}/out/{library_name}.bam.bai"
        yield from self.parent.reference_inputs().items()
        yield "targets_bed", self.parent.library_kit_paths("kits", library_name).targets

    def _get_output_files_run(self):
        prefix = "work/{library_name}/report/alfred_qc/{library_name}.alfred.json.gz"
        return {"json": prefix, "json_md5": prefix + ".md5"}

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.target_cov_report.log"}

    default_resource_usage = ResourceUsage(threads=2, runtime="4h", mem="20GB")


class BamCollectDocStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Generate depth of coverage files."""

    #: Step name
    name = "bam_collect_doc"

    #: Run for "dna" and "dna_long" only
    tool_categories = ("dna", "dna_long")

    #: Class available actions
    actions = ("run",)

    def skip_result_files_for_library(self, library_name: str) -> bool:
        return not self.config.bam_collect_doc.enabled or super().skip_result_files_for_library(
            library_name
        )

    def __init__(self, parent):
        super().__init__(parent)

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "bam", "work/{library_name}/out/{library_name}.bam".format(**wildcards)
        yield "bai", "work/{library_name}/out/{library_name}.bam.bai".format(**wildcards)
        yield from self.parent.reference_inputs().items()

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "window_length": self.config.bam_collect_doc.window_length,
        }

    @dictify
    def _get_output_files_run(self):
        prefix = "work/{library_name}/report/cov/{library_name}"
        yield "vcf", prefix + ".cov.vcf.gz"
        yield "vcf_md5", prefix + ".cov.vcf.gz.md5"
        yield "vcf_tbi", prefix + ".cov.vcf.gz.tbi"
        yield "vcf_tbi_md5", prefix + ".cov.vcf.gz.tbi.md5"
        yield "cov_bw", prefix + ".cov.bw"
        yield "cov_bw_md5", prefix + ".cov.bw.md5"
        yield "mq_bw", prefix + ".mq.bw"
        yield "mq_bw_md5", prefix + ".mq.bw.md5"

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.bam_collect_doc.log"}

    default_resource_usage = ResourceUsage(threads=1, runtime="24h", mem="2GB")


class NgsChewStepPart(ReportGetResultFilesMixin, BaseStepPart):
    """Analyze BAM File with ``ngs-chew``, e.g., ``fingerprint``"""

    #: Step name
    name = "ngs_chew"

    #: Class available actions
    actions = ("fingerprint",)

    def __init__(self, parent):
        super().__init__(parent)

    def skip_result_files_for_library(self, library_name: str) -> bool:
        return not self.config.ngs_chew_fingerprint.enabled

    @dictify
    def _get_input_files_fingerprint(self, wildcards):
        yield "bam", "work/{library_name}/out/{library_name}.bam".format(**wildcards)
        yield from self.parent.reference_inputs().items()

    def _get_output_files_fingerprint(self):
        prefix = "work/{library_name}/report/fingerprint/{library_name}.npz"
        return {"npz": prefix, "npz_md5": prefix + ".md5"}

    def _get_log_file_fingerprint(self):
        return {"log": "work/{library_name}/log/{library_name}.ngs_chew_fingerprint.log"}

    default_resource_usage = ResourceUsage(threads=1, runtime="4h", mem="2GB")


class NgsMappingWorkflow(BaseStep):
    """Perform NGS Mapping"""

    #: Step name
    name = "ngs_mapping"

    #: All libraries, background included; the tool maps those of its category (DNA, long-read
    #: DNA or RNA) and ``library_selection`` narrows them
    library_default = ""
    background_results = True

    produces = [
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
    ]

    config_model_class = NgsMappingConfigModel

    @classmethod
    def task_produces(cls, config, upstream):
        """DNA or RNA alignments, depending on the mapping tool."""
        molecule = "dna" if Tool(config.tool).is_dna() else "rna"
        return (DataSignature(DataType.ALIGNMENTS, frozenset({molecule})),)

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        """Return local BAM/BAI output paths for an alignment signature.

        Arguments:
            signature: Expected to satisfy ``DataSignature(DataType.ALIGNMENTS, …)``.
            **kwargs: Accepts ``library_name`` for concrete path rendering; falls back to
                the ``{library_name}`` wildcard placeholder.
        """
        lib = kwargs.get("library_name", "{library_name}")
        paths = {
            "bam": f"work/{lib}/out/{lib}.bam",
            "bai": f"work/{lib}/out/{lib}.bam.bai",
            "bamstats": f"work/{lib}/report/bam_qc/{lib}.bam.bamstats.txt",
            "flagstats": f"work/{lib}/report/bam_qc/{lib}.bam.flagstats.txt",
            "idxstats": f"work/{lib}/report/bam_qc/{lib}.bam.idxstats.txt",
        }
        if config.target_coverage_report and config.target_coverage_report.enabled:
            paths["alfred_qc"] = f"work/{lib}/report/alfred_qc/{lib}.alfred.json.gz"
        if config.tool == "star":
            # Strandedness inference and the gene counts it passes on
            paths["gene_counts"] = f"work/{lib}/out/{lib}.GeneCounts.tab"
            paths["strandedness"] = f"work/{lib}/strandedness/{lib}.decision.json"
        return paths

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        match self.config.tool:
            case "bwa":
                selected_mapper = BwaStepPart
            case "bwa_mem2":
                selected_mapper = BwaMem2StepPart
            case "mbcs":
                selected_mapper = MBCsStepPart
            case "minimap2":
                selected_mapper = Minimap2StepPart
            case "star":
                selected_mapper = StarStepPart
            case _:
                raise NotImplementedError(f"Unknown tool: {self.config.tool}")
        self.register_sub_step_classes(
            (
                selected_mapper,
                StrandednessStepPart,
                TargetCovReportStepPart,
                BamCollectDocStepPart,
                NgsChewStepPart,
            )
        )
        self._check_kits()
        # Validate project
        self.validate_project(config=self.config)

    @cached_property
    def library_infos(self) -> dict[str, dict[str, Any]]:
        """Return library name -> row of the library dataframe, background data sets included.

        The task maps the libraries of the background data sets as well. A name that several data
        sets have must be settled by the task's ``library_selection``.
        """
        return {name: self.library_row(name) for name in self._library_rows}

    def library_tool_category(self, library_name: str) -> str:
        """Return the tool category of a library: its extraction type, ``_long`` for long reads."""
        info = self.library_infos[library_name]
        suffix = "_long" if info["seq_platform"] in ("ONT", "PacBio") else ""
        return f"{info['extraction_type']}{suffix}"

    def get_index_path(self, tool_name: str) -> str:
        """Return the index of ``tool_name`` from the ``depends_on.index`` task, tagged with it."""
        signature = DataSignature(DataType.INDEX, frozenset({tool_name}))
        return self.get_upstream_paths("index", signature=signature).index

    def get_index_files(self, tool_name: str) -> list[str]:
        """Return the files of the mapper index, which are rule inputs.

        The ``depends_on.index`` task provides them. BWA and BWA-MEM2 indices are file prefixes with
        ``.amb`` first. The STAR index is a directory, given by its first file: the wrappers
        take its ``dirname``.
        """
        path = self.get_index_path(tool_name)
        match tool_name:
            case "bwa":
                return [path + ext for ext in BWA_INDEX_EXTENSIONS]
            case "bwa_mem2":
                return [path + ext for ext in BWA_MEM2_INDEX_EXTENSIONS]
            case "star":
                return [os.path.join(path, name) for name in STAR_INDEX_FILES]
        return [path]

    def _check_kits(self):
        """Check that ``depends_on.kits`` has the kit of each library that needs one."""
        if not self.config.depends_on.kits:
            return
        by_category = {
            name: self.library_tool_category(name) for name in self.libraries["library_name"]
        }
        if self.config.target_coverage_report and self.config.target_coverage_report.enabled:
            self.check_library_kits(
                "kits",
                [name for name, cat in by_category.items() if cat in ("dna", "dna_long")],
                require=False,
            )
        if self.config.tool == "mbcs" and self.config.mbcs.use_barcodes:
            self.check_library_kits(
                "kits", [name for name, cat in by_category.items() if cat == "dna"]
            )

    def validate_project(self, config):
        """Validates project.

        Method compares sample information included in the sample sheet and the configuration. If
        sheet contains 'DNA' samples, the configured tool must support DNA mapping. Similarly, if
        it contains 'RNA', the tool must support RNA mapping.

        :param config: SnappyStepModel with configurations as found in the project's yaml file.
        :type config: SnappyStepModel
        """
        df = self.build_library_dataframe()
        if df.empty:
            return

        # library_selection, or the default of each library's sheet kind
        df = select_libraries(df, getattr(self.config, "library_selection", None))

        if df.empty:
            return

        # Check if selected libraries contain RNA for non-STAR tools
        has_rna = (df["extraction_type"] == "rna").any()
        has_dna = (df["extraction_type"] == "dna").any()
        if has_rna and self.config.tool != "star":
            if has_dna:
                self.logger.warning(
                    "Sample sheet contains both DNA and RNA but the configured tool '%s' only "
                    "supports DNA mapping. RNA libraries will be skipped.",
                    self.config.tool,
                )
            else:
                raise InvalidConfiguration(
                    "Sample sheet contains RNA but the configured tool '{}' does not support "
                    "RNA mapping. Use 'star' for RNA mapping.".format(self.config.tool)
                )

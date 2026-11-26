"""UMI consensus mapping of DNA libraries with fgumi

The step maps the reads of libraries with unique molecular identifiers (UMIs) and calls one consensus read per molecule.
Its result files are sorted, indexed BAM files of the filtered consensus reads, which downstream tasks read like the BAM files of a mapping task.
The chain is fgumi's best-practice pipeline: ``extract`` (and ``correct``), ``fastq | bwa-mem3 mem | zipper | sort | group``, then ``duplex`` or ``simplex``, ``realign``, ``filter``, (``clip``) and the coordinate sort.
Two jobs per library run it: the grouping job streams from the FASTQ files to the grouped BAM file, the consensus job from the grouped BAM file to the final one.

The ``consensus`` field selects the library design:

- ``duplex``: UMIs on both reads (for example ``8M+T 8M+T``), the reads of both strands of a molecule form one duplex consensus read.
  Grouping uses ``--strategy paired``, the caller is ``fgumi duplex``, and the filter takes up to three ``min_reads`` values (duplex read, AB strand, BA strand) and ``require_single_strand_agreement``.
- ``simplex``: one UMI per molecule (for example ``8M+T +T``), one consensus read per UMI family.
  Grouping uses ``--strategy adjacency``, the caller is ``fgumi simplex``, and the filter takes one ``min_reads`` value.

The section ``fgumi.duplex`` or ``fgumi.simplex`` holds the settings of the caller and the mode-specific settings of the filter.
Both run with the defaults of fgumi's best practices (the caller keeps every molecule with ``min_reads`` 1 and masks bases below quality 20; the filter requires 3 raw reads for simplex and 3, 2, 2 for duplex, masks bases below quality 40 and drops reads with a raw-read error rate above 0.025 or more than 20 percent of masked bases).
Recommended settings by application, from fgumi's best practices: high sensitivity with simplex, ``filter_min_reads: 2``, ``max_base_error_rate: 0.2`` and ``max_no_call_fraction: 0.3``; high specificity with duplex, ``filter_min_reads: [10, 5, 3]`` and ``max_no_call_fraction: 0.1``; liquid biopsy with duplex, ``filter_min_reads: [3, 2, 2]`` and ``max_base_error_rate: 0.05``.

The task maps all DNA short-read libraries of its sample sheets unless ``library_selection`` narrows them, including those of background data sets.
Every library needs a read structure (see below).

Reads come from ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
The bwa-mem3 index comes from ``depends_on.index``, a ``reference_index`` task with ``tool: bwa_mem3`` (:ref:`step_reference_index`) built from the genome in ``depends_on.reference``, or an ``external_data`` task.
With a ``correct`` section, ``depends_on.umis`` names an ``external_data`` task with the kit's fixed UMIs (``produces: {type: tabular, tags: [umis]}``, ``files: {umis: <one UMI per line>}``), and ``fgumi correct`` corrects the extracted UMIs against them.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.bam`` and ``{library_name}/out/{library_name}.bam.bai``, the filtered consensus reads
- ``{library_name}/report/bam_qc/{library_name}.bam.bamstats.txt``, ``.bam.flagstats.txt`` and ``.bam.idxstats.txt``, the ``samtools`` statistics of the consensus reads
- ``{library_name}/report/umi_qc/{library_name}.group.family_sizes.txt``, ``.group.grouping_metrics.txt`` and ``.group.position_group_sizes.txt``, the metrics of ``fgumi group``
- ``{library_name}/report/umi_qc/{library_name}.consensus.family_sizes.txt``, ``.consensus.umi_counts.txt`` and ``.consensus.duplex_yield_metrics.txt`` or ``.consensus.simplex_yield_metrics.txt``, the metrics of the caller; duplex adds ``.consensus.duplex_family_sizes.txt``
- ``{library_name}/report/umi_qc/{library_name}.consensus.txt`` and ``.filter.txt``, the statistics of the caller and the filter; with ``correct`` also ``.correct.txt``, with ``clip`` also ``.clip.txt``
- ``qc_summary/out/qc_summary.tsv``, one row per library with the ``samtools stats`` counts and the UMI metrics (see "QC summary" in the notes)

The grouped BAM file ``{library_name}/out/{library_name}.grouped.bam`` (template-coordinate sorted, with the ``MI`` tags of ``fgumi group``) is an intermediate: it stays in ``work/`` to rerun the consensus with other settings, and ``delete_intermediates: [intermediates]`` removes it.

Downstream tasks read the keys ``bam``, ``bai``, ``bamstats``, ``flagstats`` and ``idxstats``.

Reads and sample sheet
----------------------

A data set or an ``external_data`` task locates the reads below its ``search_paths``, in a folder per library, as for :ref:`step_ngs_mapping`.
Each read group gets its own ``@RG`` in the grouped BAM file, with the name that the ``readgroup`` group captured as ``ID`` and ``PU``, the library name as ``SM`` and ``LB``, and ``platform`` as ``PL``.
The consensus reads of a library form one read group with the same ``SM``.

The custom field ``readStructure`` of the sample sheet gives the fgumi read structures of the FASTQ files of a read group, space separated, for example ``8M+T 8M+T`` for 8 bp UMIs at the start of both reads, or ``8M+T +T`` for a UMI on read 1 only.
A read structure is a list of segments, each a length (or ``+`` for the rest of the read) and an operator: ``T`` template, ``M`` UMI, ``S`` skip, ``B`` sample barcode, ``C`` cell barcode.
With ``consensus: duplex``, every read needs an ``M`` segment.
Libraries with another read structure or kit go into a task of their own through ``library_selection``.

Notes
-----

Fixed settings
~~~~~~~~~~~~~~

Downstream steps rely on these settings, so they are not configurable:
``bwa-mem3 mem`` runs with ``-Y`` (supplementary alignments soft clipped, so no base is lost) and ``-p`` (interleaved pairs), the realigner's bwa-mem3 as well;
``fgumi group`` records the window of the raw reads (``--realign-window``) for ``fgumi realign``;
the caller writes the per-base tags (``--output-per-base-tags``) that the filter and the clipper read;
``fgumi clip`` clips hard and keeps the per-base tags in step (``--auto-clip-attributes``);
the final sort writes the index.
The realigner's batch size is ``fgumi.realign.chunk_size`` (10 Mbp); fgumi's default of 150 Mbp deadlocks on consensus reads with per-base tags.

Local scratch
~~~~~~~~~~~~~

The grouping job writes the unmapped BAM file of each read group into the job's temporary directory, since ``fgumi zipper`` reads it twice, and with several read groups also their sorted BAM files before ``fgumi merge``.
It needs local scratch of about the size of the FASTQ files.

QC summary
~~~~~~~~~~

One job per task reads the ``samtools`` statistics and the fgumi metrics of all libraries.
It writes the summary table ``qc_summary/out/qc_summary.tsv`` with one row per library: the read counts and fractions of the consensus BAM file, the templates that grouping accepted and discarded, the families and duplexes at full depth, and the filter's pass rate.
It renders a datavzrd site below ``qc_summary/report/html/`` with the family size distributions, the duplex yield against the sequencing depth and the UMI counts, which the Snakemake report shows (``snappy report``).
"""

import re
from typing import Any

import pandas as pd
from snakemake.iocontainers import Wildcards

from snappy_pipeline.base import ConfigError
from snappy_pipeline.models import Intermediate
from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    DataSignature,
    ReportOutput,
    ResourceUsage,
)
from snappy_pipeline.workflows.abstract.protocol import DataType
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params
from snappy_pipeline.workflows.ngs_mapping import (
    BAM_REPORTS,
    BWA_MEM3_INDEX_EXTENSIONS,
    LONG_READS,
)

from .model import Consensus, Tool
from .model import UmiMapping as UmiMappingConfigModel

#: The ``fgumi group`` strategy of each consensus mode
GROUP_STRATEGY = {Consensus.duplex: "paired", Consensus.simplex: "adjacency"}

#: A read structure: segments of a length, or ``+`` for the rest of the read, and an operator
READ_STRUCTURE = re.compile(r"^(\d+[TMBCS])*(\d+|\+)[TMBCS]$")

#: The files of ``fgumi group --metrics``
GROUP_METRICS = ("family_sizes", "grouping_metrics", "position_group_sizes")

#: The inline metrics files of the consensus caller, by mode
CONSENSUS_METRICS = {
    Consensus.duplex: ("family_sizes", "duplex_family_sizes", "umi_counts", "duplex_yield_metrics"),
    Consensus.simplex: ("family_sizes", "umi_counts", "simplex_yield_metrics"),
}

#: The BAM files of a library, and the prefix of its UMI metrics
BAM = "work/{library_name}/out/{library_name}.bam"
GROUPED_BAM = "work/{library_name}/out/{library_name}.grouped.bam"
UMI_QC = "work/{library_name}/report/umi_qc/{library_name}"


def parse_read_structures(library_name: str, value) -> list[str]:
    """The read structures of ``library_name`` from the value of its ``read_structure`` column,
    one per FASTQ file of a read group; a missing or malformed value raises ``ConfigError``."""
    if value is None or pd.isna(value) or str(value).strip() in ("", "."):
        raise ConfigError(
            f"library {library_name!r} has no readStructure in its sample sheet; "
            "umi_mapping needs one, for example 8M+T 8M+T"
        )
    structures = str(value).split()
    for structure in structures:
        if not READ_STRUCTURE.match(structure):
            raise ConfigError(
                f"library {library_name!r} has the read structure {structure!r}, which is no "
                "fgumi read structure (segments of a length or + and one of T, M, S, B, C)"
            )
    return structures


def _with_md5(key: str, path: str):
    yield key, path
    yield key + "_md5", path + ".md5"


class FgumiStepPart(BaseStepPart):
    """The fgumi chain of a library in two jobs: ``grouping``, from the FASTQ files to the grouped
    BAM file, and ``consensus``, from there to the filtered, coordinate-sorted consensus reads"""

    name = "fgumi"
    actions = ("grouping", "consensus")
    wrappers = {"grouping": "fgumi/grouping", "consensus": "fgumi/consensus"}

    #: The tool's ``extra_arguments`` go to the consensus caller
    extra_arguments_action = "consensus"

    @property
    def result_keys(self) -> dict[str, tuple[str, ...]]:
        """Every output but the grouped BAM is a result file."""
        return {
            action: tuple(key for key in self.get_output_files(action) if key != "grouped_bam")
            for action in self.actions
        }

    @property
    def report_outputs(self) -> dict[str, dict[str, ReportOutput]]:
        """The ``samtools`` reports and the fgumi metrics, everything below ``report/``."""
        result = {}
        for action in self.actions:
            result[action] = {}
            for key, path in self.get_output_files(action).items():
                if key.endswith("_md5") or "/report/" not in str(path):
                    continue
                if "/bam_qc/" in str(path):
                    result[action][key] = ReportOutput("BAM QC", caption="bam_qc.rst")
                else:
                    result[action][key] = ReportOutput("UMI QC", caption="umi_qc.rst")
        return result

    #: The grouped BAM, kept to rerun the consensus with other settings
    intermediate_outputs = {"grouping": ("grouped_bam",)}

    #: The consensus BAM, which a BQSR task may process further
    bam_outputs = {"consensus": ("bam", "bam_bai")}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.fgumi.aligner.num_threads
        return ResourceUsage(threads=threads, runtime="3d", mem=f"{int(4.5 * 1024 * threads)}MB")

    def _get_input_files_grouping(self, wildcards: Wildcards) -> dict[str, Any]:
        inputs = reads_input_files(self.parent, wildcards.library_name)
        inputs["index"] = self.parent.index_files()
        inputs.update(self.parent.reference_inputs())
        if self.config.fgumi.correct is not None:
            inputs["umis"] = self.parent.get_upstream_paths("umis").umis
        return inputs

    @dictify
    def _get_output_files_grouping(self):
        yield "grouped_bam", GROUPED_BAM
        for name in GROUP_METRICS:
            yield from _with_md5(f"group_{name}", f"{UMI_QC}.group.{name}.txt")
        if self.config.fgumi.correct is not None:
            yield from _with_md5("correct_metrics", f"{UMI_QC}.correct.txt")

    def _get_log_file_grouping(self):
        return {"log": "work/{library_name}/log/{library_name}.grouping.log"}

    def _get_params_grouping(self, wildcards: Wildcards) -> dict[str, Any]:
        library = wildcards.library_name
        groups = self.parent.read_groups(library)
        structures = self.parent.read_structures(library)
        mates = 2 if groups and groups[0].right else 1
        if len(structures) != mates:
            raise ConfigError(
                f"library {library!r} has {len(structures)} read structures "
                f"({' '.join(structures)}) but {mates} FASTQ files per read group; "
                "give one read structure per file in the sheet's readStructure field"
            )
        fgumi = self.config.fgumi
        return {
            "sample_name": library,
            "read_groups": [group.name for group in groups],
            "read_structures": structures,
            "strategy": GROUP_STRATEGY[self.config.consensus],
            "num_threads": fgumi.num_threads,
            "platform": fgumi.platform,
            "aligner": fgumi.aligner.model_dump(),
            "correct": None if fgumi.correct is None else fgumi.correct.model_dump(),
            "group": fgumi.group.model_dump(),
        } | reads_params(self.parent, library)

    def _get_input_files_consensus(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "bam": GROUPED_BAM.format(library_name=wildcards.library_name),
            "index": self.parent.index_files(),
            **self.parent.reference_inputs(),
        }

    @dictify
    def _get_output_files_consensus(self):
        exts = (".bam", ".bam.bai")
        if Intermediate.bam not in self.config.delete_intermediates:
            exts += (".bam.md5", ".bam.bai.md5")  # only result files have an .md5
        for ext in exts:
            yield ext[1:].replace(".", "_"), BAM.removesuffix(".bam") + ext
        for report in BAM_REPORTS:
            yield from _with_md5(
                f"report_{report}_txt",
                f"work/{{library_name}}/report/bam_qc/{{library_name}}.bam.{report}.txt",
            )
        for name in CONSENSUS_METRICS[self.config.consensus]:
            yield from _with_md5(f"consensus_{name}", f"{UMI_QC}.consensus.{name}.txt")
        yield from _with_md5("consensus_stats", f"{UMI_QC}.consensus.txt")
        yield from _with_md5("filter_stats", f"{UMI_QC}.filter.txt")
        if self.config.fgumi.clip is not None:
            yield from _with_md5("clip_metrics", f"{UMI_QC}.clip.txt")

    def _get_log_file_consensus(self):
        return {"log": "work/{library_name}/log/{library_name}.consensus.log"}

    def _get_params_consensus(self, wildcards: Wildcards) -> dict[str, Any]:
        fgumi = self.config.fgumi
        mode = str(self.config.consensus)
        return {
            "consensus": mode,
            "num_threads": fgumi.num_threads,
            "compression_level": fgumi.compression_level,
            "aligner": fgumi.aligner.model_dump(),
            "caller": getattr(fgumi, mode).model_dump(),
            "filter": fgumi.filter.model_dump(),
            "realign": fgumi.realign.model_dump(),
            "clip": None if fgumi.clip is None else fgumi.clip.model_dump(),
        }


class UmiQcSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the ``samtools`` statistics and the fgumi metrics of all libraries"""

    report_outputs = {
        "run": {
            "html": ReportOutput("UMI QC", caption="umi_qc_summary.rst", htmlindex="index.html")
        }
    }

    source = ("fgumi", "consensus")

    @property
    def source_keys(self) -> dict[str, str | tuple[str, str]]:
        """The ``samtools`` reports and the caller's files of the consensus action, the grouping
        metrics of the grouping action; the metrics of ``correct`` and ``clip`` when they run."""
        fgumi = self.config.fgumi
        keys: dict[str, str | tuple[str, str]] = {
            report: f"report_{report}_txt" for report in BAM_REPORTS
        }
        keys |= {
            "grouping_metrics": ("grouping", "group_grouping_metrics"),
            "group_family_sizes": ("grouping", "group_family_sizes"),
            "position_group_sizes": ("grouping", "group_position_group_sizes"),
            "family_sizes": "consensus_family_sizes",
            "umi_counts": "consensus_umi_counts",
            "yield_metrics": f"consensus_{self.config.consensus}_yield_metrics",
            "consensus_stats": "consensus_stats",
            "filter_stats": "filter_stats",
        }
        if self.config.consensus == Consensus.duplex:
            keys["duplex_family_sizes"] = "consensus_duplex_family_sizes"
        if fgumi.correct is not None:
            keys["correct_metrics"] = ("grouping", "correct_metrics")
        if fgumi.clip is not None:
            keys["clip_metrics"] = "clip_metrics"
        return keys

    format = "umi"

    def summary_options(self) -> dict[str, Any]:
        return {"consensus": str(self.config.consensus)}


class UmiMappingWorkflow(BaseStep):
    """UMI consensus mapping with fgumi"""

    #: Step name
    name = "umi_mapping"

    produces = [
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna", "consensus", str(mode)}))
        for mode in Consensus
    ]

    config_model_class = UmiMappingConfigModel
    tool_parts = {Tool.fgumi: FgumiStepPart}

    #: The short-read DNA libraries
    library_requirement = f"extraction_type == 'dna' and not {LONG_READS}"

    #: All libraries, background included; ``library_selection`` narrows them
    library_default = ""

    background_results = True

    @classmethod
    def task_produces(cls, config, upstream):
        """Consensus DNA alignments, tagged with the consensus mode."""
        return (
            DataSignature(
                DataType.ALIGNMENTS, frozenset({"dna", "consensus", str(config.consensus)})
            ),
        )

    @classmethod
    def get_output_paths(cls, config, signature=None, **kwargs) -> dict[str, str]:
        lib = kwargs.get("library_name", "{library_name}")
        return {
            "bam": f"work/{lib}/out/{lib}.bam",
            "bai": f"work/{lib}/out/{lib}.bam.bai",
            "bamstats": f"work/{lib}/report/bam_qc/{lib}.bam.bamstats.txt",
            "flagstats": f"work/{lib}/report/bam_qc/{lib}.bam.flagstats.txt",
            "idxstats": f"work/{lib}/report/bam_qc/{lib}.bam.idxstats.txt",
        }

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        self.register_sub_step_classes((UmiQcSummaryStepPart,))

    def check_config(self):
        """Every library has a read structure that fits the consensus mode."""
        for library in self.libraries["library_name"]:
            structures = self.read_structures(library)
            with_umi = ["M" in structure for structure in structures]
            if self.config.consensus == Consensus.duplex and not all(with_umi):
                raise ConfigError(
                    f"library {library!r} has the read structures {' '.join(structures)}, but "
                    "consensus: duplex needs a UMI segment (M) in every read"
                )
            if not any(with_umi):
                raise ConfigError(
                    f"library {library!r} has the read structures {' '.join(structures)} without "
                    "a UMI segment (M)"
                )

    def index_files(self) -> list[str]:
        """The files of the bwa-mem3 index of ``depends_on.index``, which are rule inputs."""
        prefix = self.get_upstream_paths("index").index
        return [prefix + ext for ext in BWA_MEM3_INDEX_EXTENSIONS]

    def read_structures(self, library_name: str) -> list[str]:
        """The fgumi read structures of the FASTQ files of a read group of ``library_name``, from
        the sheet's ``readStructure`` field."""
        row = self.library_row(library_name)
        if "read_structure" not in row:
            raise ConfigError(
                f"the sample sheet of library {library_name!r} has no readStructure field; "
                "umi_mapping needs one, for example 8M+T 8M+T"
            )
        return parse_read_structures(library_name, row["read_structure"])

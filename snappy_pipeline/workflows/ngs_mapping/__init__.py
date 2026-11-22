"""Read mapping of DNA and RNA libraries

The step aligns the reads of each library to the reference genome.
Its result files are sorted, indexed BAM files.
It also computes alignment statistics.
The fields ``mask_duplicates`` and ``trim_adapters`` of the tool section turn duplicate marking and adapter trimming on or off.
``minimap2`` does neither.

The ``tool`` field selects the aligner:

- ``bwa``, ``bwa_mem2`` and ``bwa_mem3``: short-read DNA.
- ``minimap2``: long-read DNA (PacBio, Oxford Nanopore).
- ``star``: short-read RNA.
  Its result files also include gene counts and the strandedness decision.

The task maps all libraries of its sample sheets unless ``library_selection`` narrows them.
This includes the libraries of background data sets.
A task with a DNA tool maps the DNA libraries, a task with ``star`` the RNA libraries.

Reads come from ``depends_on.reads``, which names ``data_sets``, an ``external_data`` task, or an ``adapter_trimming`` task (:ref:`step_adapter_trimming`).
The mapper index comes from ``depends_on.index``, for example an ``external_data`` task or a ``reference_index`` task (:ref:`step_reference_index`).
A ``reference_index`` task must build it from the genome in ``depends_on.reference``.

Outputs
-------

The task links the following result files below ``tasks/<task>/output/``, for each library with name ``{library_name}``.
Each file has an ``.md5`` file.

- ``{library_name}/out/{library_name}.bam`` and ``{library_name}/out/{library_name}.bam.bai``
- ``{library_name}/report/bam_qc/{library_name}.bam.bamstats.txt``, ``.bam.flagstats.txt`` and ``.bam.idxstats.txt``
- ``{library_name}/out/{library_name}.cov.vcf.gz`` (with ``.tbi``), ``.cov.bw`` and ``.mq.bw``, if ``bam_collect_doc`` is enabled
- ``qc_summary/out/qc_summary.tsv``, one row per library with the ``samtools stats`` counts and fractions (see "QC summary" in the notes)
- with ``star``, for RNA libraries: ``{library_name}/out/{library_name}.GeneCounts.tab``, ``{library_name}/out/{library_name}.Junctions.tab``, ``{library_name}/out/{library_name}.Log.final.out`` (STAR's mapping summary) and ``{library_name}/strandedness/{library_name}.decision.json``

Downstream tasks read the keys ``bam``, ``bai``, ``bamstats``, ``flagstats``, and ``idxstats``, and for RNA libraries ``gene_counts`` and ``strandedness``.

Reads and sample sheet
----------------------

A data set or an ``external_data`` task locates the reads below its ``search_paths``, in a folder per library.
Each entry of ``search_patterns`` holds regular expressions for the path below that folder.
The keys ``left`` and ``right`` select the mates, and a ``(?P<readgroup>...)`` group pairs the mates of one read group.
``bwa``, ``bwa_mem2`` and ``bwa_mem3`` map each read group with its own ``@RG``.
Its ``ID`` is the name that the ``readgroup`` group captured, for example ``<library>_L001``.
Its ``SM`` is the library name, and its ``PL`` is the mapper's ``platform``.
The three BWA mappers mask duplicates (``mask_duplicates``) across all read groups of a library.

By default, every pattern must match as many ``left`` as ``right`` files.
With ``mixed_se_pe: true`` in the data set or the ``external_data`` task, ``right`` may match nothing.
Each library is then either single-end or paired-end, and libraries may differ.

Notes
-----

QC summary
~~~~~~~~~~

One job per task reads the ``samtools`` statistics of all mapped libraries.
It writes the summary table ``qc_summary/out/qc_summary.tsv`` with one row per library.
Each row holds the number of reads, the fractions of mapped, properly paired, duplicated and MAPQ 0 reads, the error rate, the read length, the base quality and the insert size.
The job also writes a site of tables and charts below ``work/qc_summary/report/html/``, which ``snappy report`` shows.
The site holds the summary, the insert size and coverage distributions, the mean base quality by cycle and the reads per contig.

STAR
~~~~

STAR needs a gene model for the gene counts and for the transcriptome alignment.
The gene model comes from ``depends_on.features``.
STAR writes counts for unstranded, forward and reverse protocols.
The step keeps the column of the library's strand and writes it as the result file ``GeneCounts.tab``.
The strand is the protocol set in ``strandedness.strand`` (0 unstranded, 1 forward, 2 reverse).
With the default ``-1``, the step infers the strand with ``infer_experiment`` from RSeQC.
``infer_experiment`` runs in both cases and reads the BED file of ``depends_on.exons``.
The file holds non-overlapping, sorted intervals of ubiquitously expressed genes with a well-defined strand, in 12 columns (6 also work).

With ``star.transcriptome: true``, STAR also writes the reads mapped to the transcriptome to a BAM file.
This BAM file has no index, because its reads are not sorted by position.

BQSR
~~~~

``alignment_processing`` with ``tool: gatk4_bqsr`` recalibrates the base qualities of the mapping.
:ref:`somatic_ngs` shows the complete chain of tasks, from adapter trimming to BQSR.
"""

import os
from typing import Any

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
from snappy_pipeline.library_table import select_libraries
from snappy_pipeline.workflows.abstract.protocol import DataType
from snappy_pipeline.workflows.common.qc_summary import QcSummaryStepPart
from snappy_pipeline.workflows.common.reads import reads_input_files, reads_params
from snappy_pipeline.workflows.common.strandedness import strandedness_params

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

from .model import NgsMapping as NgsMappingConfigModel
from .model import Tool

#: Files of a BWA index below its prefix. The wrapper takes the prefix from the first.
BWA_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt", ".pac", ".sa")

#: Files of a BWA-MEM2 index below its prefix. The wrapper takes the prefix from the first.
BWA_MEM2_INDEX_EXTENSIONS = (".amb", ".ann", ".0123", ".bwt.2bit.64", ".pac")

#: Files of a BWA-MEM3 index below its prefix (no ``.0123``). The wrapper takes the prefix from
#: the first.
BWA_MEM3_INDEX_EXTENSIONS = (".amb", ".ann", ".bwt.2bit.64", ".pac")

#: Files of a STAR index directory that mark it as the index
STAR_INDEX_FILES = ("Genome", "SA", "SAindex")

#: The BAM QC reports of ``samtools``
BAM_REPORTS = ("bamstats", "flagstats", "idxstats")

#: The long-read libraries, by the ``seqPlatform`` of the sample sheet, as a pandas query
LONG_READS = "seq_platform in ['ONT', 'PacBio']"


class ReadMappingStepPart(BaseStepPart):
    """Base class for read mapping step parts"""

    #: Class available actions. Read mapping step parts only support action "run".
    actions = ("run",)

    result_keys = ("bam", "bam_bai", "bam_md5", "bam_bai_md5") + tuple(
        f"report_{report}_txt{suffix}" for report in BAM_REPORTS for suffix in ("", "_md5")
    )

    report_outputs = {
        "run": {
            f"report_{report}_txt": ReportOutput("BAM QC", caption="bam_qc.rst")
            for report in BAM_REPORTS
        }
    }

    #: The BAM, which a BQSR task may process further
    bam_outputs = {"run": ("bam", "bam_bai")}

    #: The libraries that the tool maps, as a pandas query
    library_requirement: str

    base_path_out = "work/{{library_name}}/out/{{library_name}}{ext}"

    def _get_input_files_run(self, wildcards):
        return reads_input_files(self.parent, wildcards.library_name)

    @dictify
    def _get_output_files_run(self):
        """Return the output files that all read mapping sub steps write."""
        exts = (".bam", ".bam.bai")
        if Intermediate.bam not in self.config.delete_intermediates:
            exts += (".bam.md5", ".bam.bai.md5")  # only result files have an .md5
        for ext in exts:
            yield ext[1:].replace(".", "_"), self.base_path_out.format(ext=ext)
        for suffix in ("", ".md5"):
            for report in BAM_REPORTS:
                key = f"report_{report}_txt" + suffix.replace(".", "_")
                yield (
                    key,
                    f"work/{{library_name}}/report/bam_qc/{{library_name}}.bam.{report}.txt{suffix}",
                )

    def _get_log_file_run(self):
        return {"log": "work/{library_name}/log/{library_name}.mapping.log"}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "sample_name": wildcards.library_name,
        } | reads_params(self.parent, wildcards.library_name)

    def _read_group_names(self, wildcards: Wildcards) -> list[str]:
        """The read group IDs of the BWA mappers: the names of the discovered read groups, in
        the order of their files"""
        return [group.name for group in self.parent.read_groups(wildcards.library_name)]


class BwaStepPart(ReadMappingStepPart):
    """Short-read DNA mapping with BWA; BWA-MEM2 and BWA-MEM3 differ only in their name"""

    #: Step name
    name = "bwa"

    wrappers = {"run": "bwa"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.tool_section().num_threads_align
        return ResourceUsage(threads=threads, runtime="3d", mem=f"{int(4.5 * 1024 * threads)}MB")

    library_requirement = f"extraction_type == 'dna' and not {LONG_READS}"

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {"index": self.parent.get_index_files(self.name)}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return (
            super()._get_params_run(wildcards)
            | self.tool_section().model_dump(by_alias=True)
            | {"read_groups": self._read_group_names(wildcards)}
        )


class BwaMem2StepPart(BwaStepPart):
    name = "bwa_mem2"
    wrappers = {"run": "bwa_mem2"}


class BwaMem3StepPart(BwaStepPart):
    name = "bwa_mem3"
    wrappers = {"run": "bwa_mem3"}


class Minimap2StepPart(ReadMappingStepPart):
    """Support for performing long-read alignment using minimap2"""

    #: Step name
    name = "minimap2"

    wrappers = {"run": "minimap2"}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.minimap2.mapping_threads
        return ResourceUsage(threads=threads, runtime="2d", mem=f"{int(3.5 * threads)}GB")

    library_requirement = f"extraction_type == 'dna' and {LONG_READS}"

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {"index": self.parent.get_index_files("minimap2")}

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        params = super()._get_params_run(wildcards)
        params |= self.config.minimap2.model_dump(by_alias=True)
        params["library_kit"] = self.parent.library_row(wildcards.library_name)["library_kit"]
        params["library_name"] = wildcards.library_name
        return params


class StarStepPart(ReadMappingStepPart):
    """Support for performing NGS alignment using STAR"""

    #: Step name
    name = "star"

    wrappers = {"run": "star"}

    #: The gene counts of STAR (all strands) are input of the strandedness step, whose result
    #: file holds the counts of the inferred strand
    result_keys = ReadMappingStepPart.result_keys + (
        "junctions",
        "junctions_md5",
        "log_final",
        "log_final_md5",
        "transcriptome",
        "transcriptome_md5",
    )

    intermediate_outputs = {"run": ("gene_counts",)}

    @property
    def default_resource_usage(self) -> ResourceUsage:
        threads = self.config.star.num_threads_align
        return ResourceUsage(threads=threads, runtime="2d", mem=f"{int(3.5 * threads)}GB")

    library_requirement = f"extraction_type == 'rna' and not {LONG_READS}"

    def _get_input_files_run(self, wildcards):
        inputs = super()._get_input_files_run(wildcards)
        return inputs | {
            "features": self.parent.get_upstream_paths("features").gtf,
            "index": self.parent.get_index_files("star"),
        }

    def _get_output_files_run(self):
        """Override base class' function to make Snakemake aware of extra files for STAR."""
        output_files = super()._get_output_files_run()
        output_files["gene_counts"] = self.base_path_out.format(ext=".ReadsPerGene.tab")
        output_files["junctions"] = self.base_path_out.format(ext=".Junctions.tab")
        output_files["junctions_md5"] = output_files["junctions"] + ".md5"
        output_files["log_final"] = self.base_path_out.format(ext=".Log.final.out")
        output_files["log_final_md5"] = output_files["log_final"] + ".md5"
        if self.config.star.transcriptome:
            output_files["transcriptome"] = self.base_path_out.format(ext=".toTranscriptome.bam")
            output_files["transcriptome_md5"] = output_files["transcriptome"] + ".md5"
        return output_files

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        parent_args = super()._get_params_run(wildcards)
        parent_args.update(self.config.star.model_dump(by_alias=True))
        return parent_args


class StrandednessStepPart(BaseStepPart):
    """The strand of each RNA library, and STAR's gene counts of that strand"""

    name = "strandedness"
    actions = ("infer", "counts")
    wrappers = {"infer": "rseqc", "counts": "strandedness/counts"}

    #: ``counts`` only selects columns of the STAR counts
    extra_arguments_action = "infer"

    result_keys = {"infer": ("decision", "decision_md5"), "counts": ("counts", "counts_md5")}

    def _get_input_files_infer(self, wildcards):
        # Must reference the mapper-specific work directory
        inputs = {"bam": "work/{library_name}/out/{library_name}.bam".format(**wildcards)}
        inputs["exon_bed"] = self.parent.get_upstream_paths("exons").bed
        return inputs

    def _get_output_files_infer(self):
        tpl = "work/{library_name}/strandedness/{library_name}"
        return {
            "tsv": tpl + ".infer.txt",
            "decision": tpl + ".decision.json",
            "decision_md5": tpl + ".decision.json.md5",
        }

    def _get_log_file_infer(self):
        return {"log": "work/{library_name}/log/{library_name}.strandedness.log"}

    def _get_params_infer(self, wildcards: Wildcards) -> dict[str, Any]:
        return strandedness_params(self.config.strandedness, wildcards.library_name)

    def _get_input_files_counts(self, wildcards):
        return {
            "counts": "work/{library_name}/out/{library_name}.ReadsPerGene.tab".format(**wildcards),
            "decision": "work/{library_name}/strandedness/{library_name}.decision.json".format(
                **wildcards
            ),
        }

    def _get_output_files_counts(self):
        return {
            "counts": "work/{library_name}/out/{library_name}.GeneCounts.tab",
            "counts_md5": "work/{library_name}/out/{library_name}.GeneCounts.tab.md5",
        }

    def _get_log_file_counts(self):
        return {"log": "work/{library_name}/log/{library_name}.gene_counts.log"}


class BamCollectDocStepPart(BaseStepPart):
    """Generate depth of coverage files."""

    #: Step name
    name = "bam_collect_doc"

    #: Class available actions
    actions = ("run",)

    wrappers = {"run": "maelstrom/bam_collect_doc"}

    result_keys = (
        "vcf",
        "vcf_md5",
        "vcf_tbi",
        "vcf_tbi_md5",
        "cov_bw",
        "cov_bw_md5",
        "mq_bw",
        "mq_bw_md5",
    )

    default_resource_usage = ResourceUsage(threads=1, runtime="24h", mem="2GB")

    @dictify
    def _get_input_files_run(self, wildcards):
        yield "bam", "work/{library_name}/out/{library_name}.bam".format(**wildcards)
        yield "bai", "work/{library_name}/out/{library_name}.bam.bai".format(**wildcards)
        yield from self.parent.reference_inputs().items()

    @dictify
    def _get_output_files_run(self):
        prefix = "work/{library_name}/out/{library_name}"
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

    def _get_params_run(self, wildcards: Wildcards) -> dict[str, Any]:
        return {
            "window_length": self.config.bam_collect_doc.window_length,
        }


class BamQcSummaryStepPart(QcSummaryStepPart):
    """Tables and charts of the ``samtools`` statistics of all mapped libraries"""

    report_outputs = {
        "run": {
            "html": ReportOutput("BAM QC", caption="bam_qc_summary.rst", htmlindex="index.html")
        }
    }

    @property
    def source(self) -> tuple[str, str]:
        return (self.config.tool, "run")

    source_keys = {
        "bamstats": "report_bamstats_txt",
        "flagstats": "report_flagstats_txt",
        "idxstats": "report_idxstats_txt",
    }

    format = "samtools"


class NgsMappingWorkflow(BaseStep):
    """Perform NGS Mapping"""

    #: Step name
    name = "ngs_mapping"

    produces = [
        DataSignature(DataType.ALIGNMENTS, frozenset({"dna"})),
        DataSignature(DataType.ALIGNMENTS, frozenset({"rna"})),
    ]

    config_model_class = NgsMappingConfigModel

    tool_parts = {
        Tool.bwa: BwaStepPart,
        Tool.bwa_mem2: BwaMem2StepPart,
        Tool.bwa_mem3: BwaMem3StepPart,
        Tool.minimap2: Minimap2StepPart,
        Tool.star: StarStepPart,
    }

    #: All libraries, background included; the tool maps those of its ``library_requirement``
    #: (short-read DNA, long-read DNA or short-read RNA) and ``library_selection`` narrows them
    library_default = ""

    background_results = True

    @classmethod
    def task_library_requirement(cls, config):
        """The libraries that the tool maps."""
        return cls.tool_parts[config.tool].library_requirement

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
        if config.tool == "star":
            # Strandedness inference and the gene counts it passes on
            paths["gene_counts"] = f"work/{lib}/out/{lib}.GeneCounts.tab"
            paths["strandedness"] = f"work/{lib}/strandedness/{lib}.decision.json"
        return paths

    def __init__(self, workflow, project, task_name):
        super().__init__(workflow, project, task_name)
        if self.config.tool == Tool.star:
            # STAR's gene counts of the strand that RSeQC infers
            parts = (StrandednessStepPart,)
        elif self.config.bam_collect_doc.enabled:
            parts = (BamCollectDocStepPart,)
        else:
            parts = ()
        self.register_sub_step_classes((*parts, BamQcSummaryStepPart))
        # Validate project
        self.validate_project(config=self.config)

    def get_index_path(self, tool_name: str) -> str:
        """Return the index of ``tool_name`` from the ``depends_on.index`` task, tagged with it."""
        signature = DataSignature(DataType.INDEX, frozenset({tool_name}))
        return self.get_upstream_paths("index", signature=signature).index

    def get_index_files(self, tool_name: str) -> list[str]:
        """Return the files of the mapper index, which are rule inputs.

        The ``depends_on.index`` task provides them. BWA, BWA-MEM2 and BWA-MEM3 indices are file
        prefixes with ``.amb`` first. The STAR index is a directory, given by its first file: the
        wrappers take its ``dirname``.
        """
        path = self.get_index_path(tool_name)
        match tool_name:
            case "bwa":
                return [path + ext for ext in BWA_INDEX_EXTENSIONS]
            case "bwa_mem2":
                return [path + ext for ext in BWA_MEM2_INDEX_EXTENSIONS]
            case "bwa_mem3":
                return [path + ext for ext in BWA_MEM3_INDEX_EXTENSIONS]
            case "star":
                return [os.path.join(path, name) for name in STAR_INDEX_FILES]
        return [path]

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
                raise ConfigError(
                    "Sample sheet contains RNA but the configured tool '{}' does not support "
                    "RNA mapping. Use 'star' for RNA mapping.".format(self.config.tool)
                )

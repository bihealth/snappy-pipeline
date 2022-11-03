# -*- coding: utf-8 -*-
"""Implementation of the gCNV COHORT mode methods - used to build models.
"""

from snakemake.io import touch

from snappy_pipeline.utils import dictify
from snappy_pipeline.workflows.gcnv.gcnv_common import GcnvStepPart


class BuildGcnvModelStepPart(GcnvStepPart):
    """Class with methods to build GATK4 gCNV models"""

    #: Class available actions
    actions = (
        "preprocess_intervals",
        "annotate_gc",
        "filter_intervals",
        "scatter_intervals",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "post_germline_calls",
    )

    def __init__(self, parent):
        super().__init__(parent)

    @staticmethod
    @dictify
    def _get_input_files_annotate_gc(wildcards):
        name_pattern = "gcnv_preprocess_intervals.{wildcards.library_kit}".format(
            wildcards=wildcards
        )
        ext = "interval_list"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @dictify
    def _get_input_files_filter_intervals(self, wildcards):
        yield from self._get_input_files_annotate_gc(wildcards).items()
        name_pattern = "gcnv_annotate_gc.{wildcards.library_kit}".format(wildcards=wildcards)
        ext = "tsv"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )
        key = "covs"
        covs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                covs.append(
                    "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                        name_pattern=name_pattern, ext="tsv"
                    )
                )
        yield key, covs

    @staticmethod
    @dictify
    def _get_input_files_scatter_intervals(wildcards):
        ext = "interval_list"
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}".format(**wildcards)
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards):
        """Yield input files for ``contig_ploidy`` rule in COHORT MODE.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        ext = "interval_list"
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                tsvs.append(
                    "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                        name_pattern=name_pattern, ext=ext
                    )
                )
        yield ext, tsvs

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        """Yield input files for ``call_cnvs`` in COHORT mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        path_pattern = (
            "work/{name_pattern}/out/{name_pattern}/temp_{{shard}}/scattered.interval_list"
        )
        name_pattern = "{mapper}.gcnv_scatter_intervals.{library_kit}"
        yield "interval_list_shard", path_pattern.format(name_pattern=name_pattern)
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                path_pattern = "{mapper}.gcnv_coverage.{library_name}".format(
                    mapper=wildcards.mapper, library_name=lib
                )
                tsvs.append(
                    "work/{name_pattern}/out/{name_pattern}.{ext}".format(
                        name_pattern=path_pattern, ext=ext
                    )
                )
        yield ext, tsvs
        ext = "ploidy"
        path_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}".format(**wildcards)
        yield ext, "work/{name_pattern}/out/{name_pattern}/.done".format(name_pattern=path_pattern)
        key = "intervals"
        path_pattern = "gcnv_annotate_gc.{library_kit}"
        yield key, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=path_pattern, ext="tsv"
        )

    @staticmethod
    @dictify
    def _get_output_files_annotate_gc():
        ext = "tsv"
        name_pattern = "gcnv_annotate_gc.{library_kit}"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @staticmethod
    @dictify
    def _get_output_files_filter_intervals():
        ext = "interval_list"
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, "work/{name_pattern}/out/{name_pattern}.{ext}".format(
            name_pattern=name_pattern, ext=ext
        )

    @staticmethod
    def _get_output_files_scatter_intervals():
        return "work/{name_pattern}/out/{name_pattern}".format(
            name_pattern="{mapper}.gcnv_scatter_intervals.{library_kit}"
        )

    @staticmethod
    @dictify
    def _get_output_files_contig_ploidy():
        """Yield dictionary with output files for ``contig_ploidy`` rule in COHORT MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, touch(
            "work/{name_pattern}/out/{name_pattern}/.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )
        )

    @staticmethod
    @dictify
    def _get_output_files_call_cnvs():
        """Yield dictionary with output files for ``call_cnvs`` rle in COHORT MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        yield ext, touch(
            "work/{name_pattern}/out/{name_pattern}/.{ext}".format(
                name_pattern=name_pattern, ext=ext
            )
        )

    @staticmethod
    @dictify
    def _get_output_files_post_germline_calls():
        name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}"
        pairs = {"ratio_tsv": ".ratio.tsv", "itv_vcf": ".interval.vcf.gz", "seg_vcf": ".vcf.gz"}
        for key, ext in pairs.items():
            yield key, touch(
                "work/{name_pattern}/out/{name_pattern}{ext}".format(
                    name_pattern=name_pattern, ext=ext
                )
            )

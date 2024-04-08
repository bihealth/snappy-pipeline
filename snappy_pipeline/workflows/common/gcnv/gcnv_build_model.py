# -*- coding: utf-8 -*-
"""Implementation of the gCNV COHORT mode methods - used to build models."""

from snakemake.io import expand, touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.gcnv.gcnv_common import GcnvCommonStepPart


class AnnotateGcMixin:
    """Mixin providing functions for ``annotate_gc``"""

    @dictify
    def _get_input_files_annotate_gc(self, wildcards):
        name_pattern = f"gcnv_preprocess_intervals.{wildcards.library_kit}"
        ext = "interval_list"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    @dictify
    def _get_output_files_annotate_gc(self):
        ext = "tsv"
        name_pattern = "gcnv_annotate_gc.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    def _get_log_file_annotate_gc(self):
        name_pattern = "gcnv_annotate_gc.{library_kit}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class FilterIntervalsMixin:
    """Mixin providing functions for ``filter_intervals``"""

    @dictify
    def _get_input_files_filter_intervals(self, wildcards):
        yield from self._get_input_files_annotate_gc(wildcards).items()
        name_pattern = f"gcnv_annotate_gc.{wildcards.library_kit}"
        ext = "tsv"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"
        key = "covs"
        covs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = f"{wildcards.mapper}.gcnv_coverage.{lib}"
                ext = "tsv"
                covs.append(f"work/{name_pattern}/out/{name_pattern}.{ext}")
        yield key, covs

    @dictify
    def _get_output_files_filter_intervals(self):
        ext = "interval_list"
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    def _get_log_file_filter_intervals(self):
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class ScatterIntervalsMixin:
    """Mixin providing functions for ``scatter_intervals``"""

    @dictify
    def _get_input_files_scatter_intervals(self, wildcards):
        ext = "interval_list"
        name_pattern = f"{wildcards.mapper}.gcnv_filter_intervals.{wildcards.library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"

    def _get_output_files_scatter_intervals(self):
        name_pattern = "{mapper}.gcnv_scatter_intervals.{library_kit}"
        return f"work/{name_pattern}/out/{name_pattern}"

    def _get_log_file_scatter_intervals(self):
        name_pattern = "{mapper}.gcnv_scatter_intervals.{library_kit}"
        return f"work/{name_pattern}/log/{name_pattern}.log"


class ContigPloidyMixin:
    """Mixin providing functions for ``contig_ploidy``"""

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards):
        """Yield input files for ``contig_ploidy`` rule in COHORT MODE.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        ext = "interval_list"
        name_pattern = "{mapper}.gcnv_filter_intervals.{library_kit}"
        yield ext, f"work/{name_pattern}/out/{name_pattern}.{ext}"
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                name_pattern = f"{wildcards.mapper}.gcnv_coverage.{lib}"
                tsvs.append(f"work/{name_pattern}/out/{name_pattern}.{ext}")
        yield ext, tsvs
        # Yield path to pedigree file
        peds = []
        for library_name in sorted(self.index_ngs_library_to_pedigree):
            name_pattern = f"write_pedigree.{library_name}"
            peds.append(f"work/{name_pattern}/out/{library_name}.ped")
        yield "ped", peds

    @dictify
    def _get_output_files_contig_ploidy(self):
        """Yield dictionary with output files for ``contig_ploidy`` rule in COHORT MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_contig_ploidy.{library_kit}"
        yield ext, touch(f"work/{name_pattern}/out/{name_pattern}/.{ext}")


class CallCnvsMixin:
    """Mixin providing functions for ``call_cnvs``"""

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        """Yield input files for ``call_cnvs`` in COHORT mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        name_pattern = "{mapper}.gcnv_scatter_intervals.{library_kit}"
        path_pattern = (
            f"work/{name_pattern}/out/{name_pattern}/temp_{{shard}}/scattered.interval_list"
        )
        yield "interval_list_shard", path_pattern
        ext = "tsv"
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                path_pattern = f"{wildcards.mapper}.gcnv_coverage.{lib}"
                tsvs.append(f"work/{path_pattern}/out/{path_pattern}.{ext}")
        yield ext, tsvs
        ext = "ploidy"
        path_pattern = f"{wildcards.mapper}.gcnv_contig_ploidy.{wildcards.library_kit}"
        yield ext, f"work/{path_pattern}/out/{path_pattern}/.done"
        key = "intervals"
        path_pattern = "gcnv_annotate_gc.{library_kit}"
        yield key, f"work/{path_pattern}/out/{path_pattern}.tsv"

    @dictify
    def _get_output_files_call_cnvs(self):
        """Yield dictionary with output files for ``call_cnvs`` rle in COHORT MODE."""
        ext = "done"
        name_pattern = "{mapper}.gcnv_call_cnvs.{library_kit}.{shard}"
        yield ext, touch(f"work/{name_pattern}/out/{name_pattern}/.{ext}")


class PostGermlineCallsMixin:
    """Mixin providing functions for ``post_germline_calls``"""

    @dictify
    def _get_output_files_post_germline_calls(self):
        name_pattern = "{mapper}.gcnv_post_germline_calls.{library_name}"
        pairs = {"ratio_tsv": ".ratio.tsv", "itv_vcf": ".interval.vcf.gz", "seg_vcf": ".vcf.gz"}
        for key, ext in pairs.items():
            yield key, touch(f"work/{name_pattern}/out/{name_pattern}{ext}")


class BuildGcnvModelStepPart(
    AnnotateGcMixin,
    FilterIntervalsMixin,
    ScatterIntervalsMixin,
    ContigPloidyMixin,
    CallCnvsMixin,
    PostGermlineCallsMixin,
    GcnvCommonStepPart,
):
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
    )

    @listify
    def get_result_files(self):
        """Return list of concrete output paths"""

        # Get list with all result path template strings.  This is done using the function generating
        # the output files for the post germline calls step (will create coverage and ploidy models).
        #
        # NB: the conversion to ``str`` is necessary here as we use ``touch()`` above.
        result_path_tpls = list(map(str, self._get_output_files_post_germline_calls().values()))

        # Generate output files for all mappers and library names.
        for _, pedigree in self.index_ngs_library_to_pedigree.items():
            library_names = [
                donor.dna_ngs_library.name for donor in pedigree.donors if donor.dna_ngs_library
            ]
            for path_tpl in result_path_tpls:
                yield from expand(
                    path_tpl,
                    mapper=self.w_config["step_config"]["ngs_mapping"]["tools"]["dna"],
                    library_name=library_names,
                )

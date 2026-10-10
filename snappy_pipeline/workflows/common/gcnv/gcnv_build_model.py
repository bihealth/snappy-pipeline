# -*- coding: utf-8 -*-
"""Implementation of the gCNV COHORT mode methods - used to build models."""

from snakemake.io import expand, touch

from snappy_pipeline.utils import dictify, listify
from snappy_pipeline.workflows.common.gcnv.gcnv_common import GcnvCommonStepPart


class AnnotateGcMixin:
    """Mixin providing functions for ``annotate_gc``"""

    @dictify
    def _get_input_files_annotate_gc(self, wildcards):
        kit = wildcards.library_kit
        yield "interval_list", f"work/{kit}/out/{kit}.interval_list"
        yield "reference", self.parent.get_upstream_paths("reference").fasta

    def _get_params_annotate_gc(self, wildcards):
        gcnv = self.config.gcnv
        params = {"path_uniquely_mapable_bed": gcnv.path_uniquely_mapable_bed}
        if hasattr(gcnv, "path_target_interval_list_mapping"):  # targeted sequencing only
            mapping = gcnv.path_target_interval_list_mapping
            params["path_target_interval_list_mapping"] = [item.model_dump() for item in mapping]
        return params

    @dictify
    def _get_output_files_annotate_gc(self):
        yield "tsv", "work/{library_kit}/out/{library_kit}.annotate_gc.tsv"

    def _get_log_file_annotate_gc(self):
        return "work/{library_kit}/log/{library_kit}.annotate_gc.log"


class FilterIntervalsMixin:
    """Mixin providing functions for ``filter_intervals``"""

    @dictify
    def _get_input_files_filter_intervals(self, wildcards):
        kit = wildcards.library_kit
        yield "interval_list", f"work/{kit}/out/{kit}.interval_list"
        yield "tsv", f"work/{kit}/out/{kit}.annotate_gc.tsv"
        covs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                covs.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "covs", covs

    @dictify
    def _get_output_files_filter_intervals(self):
        yield "interval_list", "work/{library_kit}/out/{library_kit}.filter_intervals.interval_list"

    def _get_log_file_filter_intervals(self):
        return "work/{library_kit}/log/{library_kit}.filter_intervals.log"


class ScatterIntervalsMixin:
    """Mixin providing functions for ``scatter_intervals``"""

    @dictify
    def _get_input_files_scatter_intervals(self, wildcards):
        kit = wildcards.library_kit
        yield "interval_list", f"work/{kit}/out/{kit}.filter_intervals.interval_list"

    def _get_output_files_scatter_intervals(self):
        return "work/{library_kit}/out/{library_kit}.scatter_intervals"

    def _get_log_file_scatter_intervals(self):
        return "work/{library_kit}/log/{library_kit}.scatter_intervals.log"


class ContigPloidyMixin:
    """Mixin providing functions for ``contig_ploidy``"""

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards):
        """Yield input files for ``contig_ploidy`` rule in COHORT MODE.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        kit = wildcards.library_kit
        yield "interval_list", f"work/{kit}/out/{kit}.filter_intervals.interval_list"
        if self.config.gcnv.path_par_intervals:  # PAR regions to exclude
            yield "par_intervals", self.config.gcnv.path_par_intervals
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                tsvs.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "tsv", tsvs
        # Yield path to pedigree file
        peds = []
        for library_name in sorted(self.index_ngs_library_to_pedigree):
            name_pattern = f"write_pedigree.{library_name}"
            peds.append(f"work/{name_pattern}/out/{library_name}.ped")
        yield "ped", peds

    @dictify
    def _get_output_files_contig_ploidy(self):
        """Yield dictionary with output files for ``contig_ploidy`` rule in COHORT MODE."""
        yield "done", touch("work/{library_kit}/out/{library_kit}.contig_ploidy/.done")


class CallCnvsMixin:
    """Mixin providing functions for ``call_cnvs``"""

    @dictify
    def _get_input_files_call_cnvs(self, wildcards):
        """Yield input files for ``call_cnvs`` in COHORT mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        yield (
            "interval_list_shard",
            "work/{library_kit}/out/{library_kit}.scatter_intervals/temp_{shard}"
            "/scattered.interval_list",
        )
        tsvs = []
        for lib in sorted(self.index_ngs_library_to_donor):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                tsvs.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "tsv", tsvs
        kit = wildcards.library_kit
        yield "ploidy", f"work/{kit}/out/{kit}.contig_ploidy/.done"
        yield "intervals", "work/{library_kit}/out/{library_kit}.annotate_gc.tsv"

    @dictify
    def _get_output_files_call_cnvs(self):
        """Yield dictionary with output files for ``call_cnvs`` rle in COHORT MODE."""
        yield "done", touch("work/{library_kit}/out/{library_kit}.{shard}.call_cnvs/.done")


class PostGermlineCallsMixin:
    """Mixin providing functions for ``post_germline_calls``"""

    @dictify
    def _get_output_files_post_germline_calls(self):
        prefix = "work/{library_name}/out/{library_name}.post_germline_calls"
        pairs = {"ratio_tsv": ".ratio.tsv", "itv_vcf": ".interval.vcf.gz", "seg_vcf": ".vcf.gz"}
        for key, ext in pairs.items():
            yield key, touch(f"{prefix}{ext}")


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
                    library_name=library_names,
                )

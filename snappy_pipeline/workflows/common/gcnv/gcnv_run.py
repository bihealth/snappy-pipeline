# -*- coding: utf-8 -*-
"""Implementation of the gCNV CASE mode run methods."""

import re
import warnings
from pathlib import Path

from snakemake.io import expand, touch
from snakemake.iocontainers import Wildcards

from snappy_pipeline.utils import dictify, flatten, listify
from snappy_pipeline.workflows.common.gcnv.gcnv_common import (
    GcnvCommonStepPart,
    InconsistentLibraryKitsWarning,
    build_ngs_library_to_kit,
)

#: Predefined file name keys and extensions for the overall results files.
RESULT_EXTENSIONS = {
    "vcf": ".vcf.gz",
    "vcf_tbi": ".vcf.gz.tbi",
}

#: The published results files also get md5 files.
PUBLISHED_RESULT_EXTENSIONS = {
    **RESULT_EXTENSIONS,
    "vcf_md5": ".vcf.gz.md5",
    "vcf_tbi_md5": ".vcf.gz.tbi.md5",
}


class ContigPloidyMixin:
    """Methods for ``contig_ploidy``."""

    @dictify
    def _get_input_files_contig_ploidy(self, wildcards: Wildcards):
        """Yield input files for ``contig_ploidy`` rule in CASE MODE using the model task.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        """
        kit = wildcards.library_kit
        tsvs = []
        for lib in sorted(
            self.parent.primary_dna_libraries(include_background=True)["library_name"]
        ):
            if self.ngs_library_to_kit.get(lib) == kit:
                tsvs.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "tsv", tsvs
        # Yield path to pedigree file
        peds = []
        for library_name in sorted(self.parent.pedigrees(include_background=True)):
            name_pattern = f"write_pedigree.{library_name}"
            peds.append(f"work/{name_pattern}/out/{library_name}.ped")
        yield "ped", peds
        yield "model", self.parent.get_upstream_paths("gcnv_model", library_kit=kit).ploidy

    @dictify
    def _get_output_files_contig_ploidy(self):
        """Yield dictionary with output files for ``contig_ploidy`` rule in CASE MODE."""
        yield "done", touch("work/{library_kit}/out/{library_kit}.contig_ploidy/.done")


class CallCnvsMixin:
    """Methods for ``call_cnvs``."""

    @dictify
    def _get_input_files_call_cnvs(self, wildcards: Wildcards):
        """Yield input files for ``call_cnvs`` in CASE mode.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'mapper' (e.g., 'bwa')
        and 'library_kit' (e.g., 'Agilent_SureSelect_Human_All_Exon_V6').
        :type wildcards: snakemake.io.Wildcards
        """
        # Yield coverage tsv files for all library associated with kit
        coverage_files = []
        for lib in sorted(
            self.parent.primary_dna_libraries(include_background=True)["library_name"]
        ):
            if self.ngs_library_to_kit.get(lib) == wildcards.library_kit:
                coverage_files.append(f"work/{lib}/out/{lib}.coverage.tsv")
        yield "tsv", coverage_files

        # Yield ploidy files
        kit = wildcards.library_kit
        yield "ploidy", f"work/{kit}/out/{kit}.contig_ploidy/.done"

        # Yield the model of the shard
        model = self.parent.get_upstream_paths("gcnv_model", library_kit=kit, shard=wildcards.shard)
        yield "model", model.calls

    @dictify
    def _get_output_files_call_cnvs(self):
        """Yield dictionary with output files for ``call_cnvs`` rule in CASE MODE."""
        yield "done", touch("work/{library_kit}/out/{library_kit}.{shard}.call_cnvs/.done")


class PostGermlineCallsMixin:
    """Methods for ``post_germline_calls``."""

    @dictify
    def _get_output_files_post_germline_calls(self):
        prefix = "work/{library_name}/out/{library_name}.post_germline_calls"
        extensions = {
            "ratio_tsv": ".ratio.tsv",
            "itv_vcf": ".interval.vcf.gz",
            "seg_vcf": ".vcf.gz",
        }
        for key, ext in extensions.items():
            yield key, touch(f"{prefix}{ext}")

    @dictify
    def _get_input_files_post_germline_calls(self, wildcards, checkpoints):
        """Yield input files for ``post_germline_calls`` in CASE mode.

        The shards of the model task come from its checkpoint, through the local checkpoint
        ``<step>_gcnv_shards`` (a Snakemake module cannot use the checkpoints of another one).
        The Snakefile calls this function with Snakemake's ``checkpoints`` object.

        :param wildcards: Snakemake wildcards associated with rule, namely: 'library_name'
        (e.g., 'P001-N1-DNA1-WGS1').
        :type wildcards: snakemake.io.Wildcards
        """
        library_kit = self.ngs_library_to_kit.get(wildcards.library_name)
        shards_file = getattr(checkpoints, f"{self.parent.name}_gcnv_shards").get(
            library_kit=library_kit
        )
        shards = Path(shards_file.output[0]).read_text().split()

        # Yield cnv calls output
        yield (
            "calls",
            [f"work/{library_kit}/out/{library_kit}.{shard}.call_cnvs/.done" for shard in shards],
        )

        # Yield contig-ploidy output
        yield "ploidy", f"work/{library_kit}/out/{library_kit}.contig_ploidy/.done"

        # Yield the model shards (same order as the calls)
        yield (
            "model",
            [
                self.parent.get_upstream_paths(
                    "gcnv_model", library_kit=library_kit, shard=shard
                ).calls
                for shard in shards
            ],
        )


class JointGermlineCnvSegmentationMixin:
    """Methods for ``post_germline_calls``."""

    @dictify
    def _get_output_files_joint_germline_cnv_segmentation(self):
        prefix = "work/{library_name}/out/{library_name}.{kit}.joint_segmentation"
        work_files = {}
        for key, suffix in RESULT_EXTENSIONS.items():
            work_files[key] = f"{prefix}{suffix}"
        yield from work_files.items()

    def _get_log_file_joint_germline_cnv_segmentation(self):
        """Return log file **pattern** for the step ``joint_germline_cnv_segmentation``."""
        return {"log": "work/{library_name}/log/{library_name}.{kit}.joint_segmentation.log"}

    @dictify
    def _get_input_files_joint_germline_cnv_segmentation(self, wildcards):
        # Yield list of paths to input VCF files
        ped_ngs_library_names = [
            library_name
            for library_name in self.parent.pedigrees(include_background=True)[
                wildcards.library_name
            ]
            if self.ngs_library_to_kit[library_name] == wildcards.kit
        ]
        vcfs = []
        for library_name in sorted(ped_ngs_library_names):
            vcfs.append(f"work/{library_name}/out/{library_name}.post_germline_calls.vcf.gz")
        yield "vcf", vcfs
        # Yield path to interval list file
        kit = wildcards.kit
        yield "interval_list", self._interval_list_path(kit)
        # Yield path to pedigree file
        name_pattern = f"write_pedigree.{wildcards.library_name}"
        yield "ped", f"work/{name_pattern}/out/{wildcards.library_name}.ped"
        yield from self.parent.reference_inputs().items()


class MergeMultikitFamiliesMixin:
    """Methods for merging families with multiple kits.

    In the case of a single kit per family, we can simply copy the input
    to the output.
    """

    @dictify
    def _get_output_files_merge_multikit_families(self):
        prefix = "work/{library_name}/out/{library_name}"
        work_files = {}
        for key, suffix in PUBLISHED_RESULT_EXTENSIONS.items():
            work_files[key] = f"{prefix}{suffix}"
        yield from work_files.items()

    def _get_log_file_merge_multikit_families(self):
        """Return log file **pattern** for the step ``merge_multikit_families``."""
        return {"log": "work/{library_name}/log/{library_name}.merge_multikit_families.log"}

    @dictify
    def _get_input_files_merge_multikit_families(self, wildcards):
        # Obtain all kits that we need to merge for this family.
        kits = []
        for library_name in self.parent.pedigrees(include_background=True)[wildcards.library_name]:
            kit = self.ngs_library_to_kit[library_name]
            if kit not in kits:
                kits.append(kit)
        # Yield list of paths to input VCF files
        vcfs = []
        library_name = wildcards.library_name
        for kit in kits:
            vcfs.append(f"work/{library_name}/out/{library_name}.{kit}.joint_segmentation.vcf.gz")
        yield "vcf", vcfs
        yield "vcf_tbi", [f"{vcf}.tbi" for vcf in vcfs]


class ShardsMixin:
    """Methods for ``shards``, the local checkpoint that lists the shards of the model task."""

    @dictify
    def _get_input_files_shards(self, wildcards):
        model = self.parent.get_upstream_paths("gcnv_model", library_kit=wildcards.library_kit)
        yield "shards", model.shards

    @dictify
    def _get_output_files_shards(self):
        yield "shards", "work/{library_kit}/out/{library_kit}.shards.txt"

    def write_shards(self, input, output):
        """Write the names of the shard directories ``temp_{shard}`` of the model task, sorted."""
        directories = Path(input.shards).glob("temp_*")
        shards = sorted(path.name.removeprefix("temp_") for path in directories if path.is_dir())
        Path(output.shards).write_text("".join(f"{shard}\n" for shard in shards))


class RunGcnvStepPart(
    ShardsMixin,
    ContigPloidyMixin,
    CallCnvsMixin,
    PostGermlineCallsMixin,
    JointGermlineCnvSegmentationMixin,
    MergeMultikitFamiliesMixin,
    GcnvCommonStepPart,
):
    """Class with methods to run GATK4 gCNV calling in CASE mode"""

    #: Class available actions
    actions = (
        "shards",
        "coverage",
        "contig_ploidy",
        "call_cnvs",
        "post_germline_calls",
        "joint_germline_cnv_segmentation",
        "merge_multikit_families",
    )

    temp_outputs = {"call_cnvs": ("done",)}

    intermediate_outputs = {
        "coverage": ("tsv",),
        "contig_ploidy": ("done",),
        "post_germline_calls": ("seg_vcf",),
        "joint_germline_cnv_segmentation": ("vcf", "vcf_tbi"),
    }

    def __init__(self, parent):
        super().__init__(parent)
        self.ngs_library_to_kit = build_ngs_library_to_kit(
            self.parent.primary_dna_libraries(include_background=True),
            self.parent.get_task_config("gcnv_model").intervals,
            self.parent.library_kits,
        )

    def _model_config(self):
        """Return the config of the ``gcnv_model`` task."""
        return self.parent.get_task_config("gcnv_model")

    def _interval_list_path(self, library_kit: str) -> str:
        """Return the interval list of the model task, which also fixes the intervals of the calls."""
        return self.parent.get_upstream_paths("gcnv_model", library_kit=library_kit).interval_list

    @listify
    def get_result_files(self):
        """Return list of **concrete** paths to result files for the given configuration and sample sheets.

        The function will skip pedigrees where samples have inconsistent library kits and print a warning.
        """
        if self.config.tool != "gcnv":
            return

        # Get list with all result path template strings.
        result_path_tpls = [
            re.sub(r"^work/", "output/", path)
            for path in flatten(self._get_output_files_merge_multikit_families().values())
        ]

        # Iterate over all pedigrees.  Check library kit consistency for each pedigree.  In case of inconsistent
        # library kits within one pedigree, raise a warning about the quality in this pedigree.
        for (
            index_library_name,
            library_names,
        ) in self.parent.pedigrees(include_background=True).items():
            library_kits = [self.ngs_library_to_kit[library_name] for library_name in library_names]
            if len(set(library_kits)) != 1:
                names_kits = list(zip(library_names, library_kits))
                msg = (
                    "Found inconsistent library kits (more than one kit!) for pedigree with index "
                    f"{index_library_name}.  The library name/kit pairs are {names_kits}.  SNAPPY "
                    "will attempt to merge the results but the result may contain artifacts."
                )
                warnings.warn(InconsistentLibraryKitsWarning(msg))
            # Yield all concrete output paths, replacing prefix "work/" by "output/".
            for path_tpl in result_path_tpls:
                yield from expand(
                    path_tpl,
                    library_name=[index_library_name],
                )

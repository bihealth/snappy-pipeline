"""SV calling for targeted sequencing

Based on the output of ``ngs_mapping``, call structural variants from depth of coverage,
read pair, and split read signal.
"""

from copy import deepcopy
from itertools import chain
import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_background, is_not_background
from snakemake.io import expand

from snappy_pipeline.utils import DictQuery, dictify, flatten, listify
from snappy_pipeline.workflows.abstract import (
    BaseStep,
    BaseStepPart,
    LinkInStep,
    WritePedigreeStepPart,
)
from snappy_pipeline.workflows.abstract.common import (
    ForwardResourceUsageMixin,
    ForwardSnakemakeFilesMixin,
)
from snappy_pipeline.workflows.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow
from snappy_wrappers.resource_usage import ResourceUsage

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available SV callers
SV_CALLERS = ("gcnv", "delly2", "manta")

#: Minimum number of samples per kit to apply gCNV calling criteria to be analyzed
GCNV_MIN_KIT_SAMPLES = 10

#: Default configuration for the sv_calling_targeted step
DEFAULT_CONFIG = r"""
# Default configuration sv_calling_targeted
step_config:
  sv_calling_targeted:
    # Path to the ngs_mapping step.
    path_ngs_mapping: ../ngs_mapping

    # List of used tools
    tools: [gcnv, delly2, manta]  # REQUIRED

    # The following allows to define one or more set of target intervals.  This is only used by gcnv.
    #
    # Example:
    #
    # - name: "Agilent SureSelect Human All Exon V6"
    #   pattern: "Agilent SureSelect Human All Exon V6.*"
    #   path: "path/to/targets.bed"
    path_target_interval_list_mapping: []

    gcnv:
      # Path to gCNV model - will execute analysis in CASE MODE.
      #
      # Example:
      #
      # - library: "Agilent SureSelect Human All Exon V6"  # Kit name, match in path_target_interval_list_mapping
      #   contig_ploidy: /path/to/ploidy-model         # Output from `DetermineGermlineContigPloidy`
      #   model_pattern: /path/to/model_*              # Output from `GermlineCNVCaller`
      precomputed_model_paths: []

    delly2:
      path_exclude_tsv: null  # optional
      max_threads: 16
      map_qual: 1
      geno_qual: 5
      qual_tra: 20
      mad_cutoff: 9

    manta:
      max_threads: 16
"""


class GcnvTargetedStepPart(RunGcnvStepPart):
    """Targeted seq. CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self.parent.ngs_library_to_kit


class SvCallingTargetedGetResultFilesMixin:
    """Mixin that provides ``get_result_files()`` for SV calling steps"""

    @listify
    def get_result_files(self):
        """Return list of concrete output paths in ``output/``.

        The implementation will return a list of all paths with prefix ``output/` that are
        returned by ``self.get_output_files()`` for all actions in ``self.actions``.
        """
        ngs_mapping_config = DictQuery(self.w_config).get("step_config/ngs_mapping")
        for mapper in ngs_mapping_config["tools"]["dna"]:
            # Get list of result path templates.
            output_files_tmp = self.get_output_files(self.actions[-1])
            if isinstance(output_files_tmp, dict):
                output_files = output_files_tmp.values()
            else:
                output_files = output_files_tmp
            result_paths_tpls = list(
                filter(
                    lambda p: p.startswith("output/"),
                    flatten(output_files),
                )
            )
            #: Generate all concrete output paths.
            for path_tpl in result_paths_tpls:
                for library_name in self.index_ngs_library_to_pedigree.keys():
                    yield from expand(path_tpl, mapper=[mapper], library_name=library_name)


class Delly2StepPart(
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingTargetedGetResultFilesMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Delly2"""

    name = "delly2"
    actions = ("call", "merge_calls", "genotype", "merge_genotypes")

    _cheap_resource_usage = ResourceUsage(
        threads=2,
        time="4-00:00:00",
        memory=f"{7 * 1024 * 2}M",
    )
    _normal_resource_usage = ResourceUsage(
        threads=2,
        time="7-00:00:00",  # 7 days
        memory=f"{20 * 1024 * 2}M",
    )
    resource_usage_dict = {
        "call": _normal_resource_usage,
        "merge_calls": _cheap_resource_usage,
        "genotype": _normal_resource_usage,
        "merge_genotypes": _cheap_resource_usage,
    }

    def __init__(self, parent):
        super().__init__(parent)

        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

        self.donor_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.donor_ngs_library_to_pedigree.update(sheet.donor_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_call(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        token = f"{wildcards.mapper}.{wildcards.library_name}"
        yield "bam", ngs_mapping(f"output/{token}/out/{token}.bam")

    @dictify
    def _get_output_files_call(self):
        infix = "{mapper}.delly2_call.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_calls(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.delly2_call.{donor.dna_ngs_library.name}"
                bcfs.append(f"work/{infix}/out/{infix}.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_calls(self):
        infix = "{mapper}.delly2_merge_calls.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_genotype(self, wildcards):
        yield from self._get_input_files_call(wildcards).items()
        pedigree = self.donor_ngs_library_to_pedigree[wildcards.library_name]
        infix = f"{wildcards.mapper}.delly2_merge_calls.{pedigree.index.dna_ngs_library.name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"

    @dictify
    def _get_output_files_genotype(self):
        infix = "{mapper}.delly2_genotype.{library_name}"
        yield "bcf", f"work/{infix}/out/{infix}.bcf"
        yield "bcf_md5", f"work/{infix}/out/{infix}.bcf.md5"
        yield "bcf_csi", f"work/{infix}/out/{infix}.bcf.csi"
        yield "bcf_csi_md5", f"work/{infix}/out/{infix}.bcf.csi.md5"

    @dictify
    def _get_input_files_merge_genotypes(self, wildcards):
        bcfs = []
        pedigree = self.index_ngs_library_to_pedigree[wildcards.library_name]
        for donor in pedigree.donors:
            if donor.dna_ngs_library:
                infix = f"{wildcards.mapper}.delly2_genotype.{donor.dna_ngs_library.name}"
                bcfs.append(f"work/{infix}/out/{infix}.bcf")
        yield "bcf", bcfs

    @dictify
    def _get_output_files_merge_genotypes(self):
        infix = "{mapper}.delly2.{library_name}"
        work_files = {}
        work_files["vcf"] = f"work/{infix}/out/{infix}.vcf.gz"
        work_files["vcf_md5"] = f"work/{infix}/out/{infix}.vcf.gz.md5"
        work_files["vcf_tbi"] = f"work/{infix}/out/{infix}.vcf.gz.tbi"
        work_files["vcf_tbi_md5"] = f"work/{infix}/out/{infix}.vcf.gz.tbi.md5"
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(
                work_files.values(), self.get_log_file("merge_genotypes").values()
            )
        ]

    @dictify
    def get_log_file(self, action):
        """Return dict of log files in the "log" directory"""
        _ = action
        if action != "merge_genotypes":
            infix = f"delly2_{action}"
        else:
            infix = "delly2"
        prefix = f"work/{{mapper}}.{infix}.{{library_name}}/log/{{mapper}}.{infix}.{{library_name}}.sv_calling"
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


class MantaStepPart(
    ForwardSnakemakeFilesMixin,
    ForwardResourceUsageMixin,
    SvCallingTargetedGetResultFilesMixin,
    BaseStepPart,
):
    """Perform SV calling on exomes using Manta"""

    name = "manta"
    actions = ("run",)

    resource_usage_dict = {
        "run": ResourceUsage(
            threads=16,
            time="2-00:00:00",
            memory=f"{int(3.75 * 1024 * 16)}M",
        )
    }

    def __init__(self, parent):
        super().__init__(parent)
        #: Shortcuts from index NGS library name to Pedigree
        self.index_ngs_library_to_pedigree = {}
        for sheet in self.parent.shortcut_sheets:
            self.index_ngs_library_to_pedigree.update(sheet.index_ngs_library_to_pedigree)

    @dictify
    def _get_input_files_run(self, wildcards):
        ngs_mapping = self.parent.sub_workflows["ngs_mapping"]
        bams = []
        for donor in self.index_ngs_library_to_pedigree[wildcards.library_name].donors:
            if donor.dna_ngs_library:
                token = f"{wildcards.mapper}.{donor.dna_ngs_library.name}"
                bams.append(ngs_mapping(f"output/{token}/out/{token}.bam"))
        yield "bam", bams

    @dictify
    def _get_output_files_run(self):
        work_files = {}
        for name, ext in zip(EXT_NAMES, EXT_VALUES):
            infix = "{mapper}.manta.{library_name}"
            work_files[name] = f"work/{infix}/out/{infix}{ext}"
        yield from work_files.items()
        yield "output_links", [
            re.sub(r"^work/", "output/", work_path)
            for work_path in chain(work_files.values(), self.get_log_file("run").values())
        ]

    @dictify
    def get_log_file(self, action):
        """Return dict of log files in the "log" directory"""
        _ = action
        caller = self.__class__.name
        prefix = f"work/{{mapper}}.{caller}.{{library_name}}/log/{{mapper}}.{caller}.{{library_name}}.sv_calling"
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


class SvCallingTargetedWorkflow(BaseStep):
    """Perform germline targeted sequencing CNV calling"""

    #: Workflow name
    name = "sv_calling_targeted"

    sheet_shortcut_class = GermlineCaseSheet

    def __init__(self, workflow, config, config_lookup_paths, config_paths, workdir):
        super().__init__(
            workflow,
            config,
            config_lookup_paths,
            config_paths,
            workdir,
            (NgsMappingWorkflow,),
        )
        # Build mapping from NGS library name to kit
        self.ngs_library_to_kit = self._build_ngs_library_to_kit()
        # Register sub step classes so the sub steps are available
        self.register_sub_step_classes(
            (WritePedigreeStepPart, GcnvTargetedStepPart, Delly2StepPart, MantaStepPart)
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Build dictionary with sample count per library kit
        _, _, self.library_kit_counts_dict = self.pick_kits_and_donors()

    @dictify
    def _build_ngs_library_to_kit(self):
        config = DictQuery(self.w_config).get("step_config/sv_calling_targeted")
        if not config["path_target_interval_list_mapping"]:
            # No mapping given, we will use the "default" one for all.
            for donor in self.all_donors():
                if donor.dna_ngs_library:
                    yield donor.dna_ngs_library.name, "default"

        # Build mapping
        regexes = {
            item["pattern"]: item["name"] for item in config["path_target_interval_list_mapping"]
        }
        result = {}
        for donor in self.all_donors():
            if donor.dna_ngs_library and donor.dna_ngs_library.extra_infos.get("libraryKit"):
                library_kit = donor.dna_ngs_library.extra_infos.get("libraryKit")
                for pattern, name in regexes.items():
                    if re.match(pattern, library_kit):
                        yield donor.dna_ngs_library.name, name
        return result

    @classmethod
    def default_config_yaml(cls):
        """Default configuration.

        :return: Returns default config YAML, to be overwritten by project-specific one.
        """
        return DEFAULT_CONFIG

    def get_library_count(self, library_kit):
        """Get library count.

        :param library_kit: Library kit name.
        :type library_kit: str

        :return: Returns number of samples with inputted library kit. If library name not defined,
        it returns zero.
        """
        return self.library_kit_counts_dict.get(library_kit, 0)

    @listify
    def all_donors(self, include_background=True):
        """Get all donors.

        :param include_background: Boolean flag to defined if background should be included or not.
        Default: True, i.e., background will be included.

        :return: Returns list of all donors in sample sheet.
        """
        sheets = self.shortcut_sheets
        if not include_background:
            sheets = list(filter(is_not_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def all_background_donors(self):
        """Get all background donors.

        :return: Returns list of all background donors in sample sheets.
        """
        sheets = deepcopy(self.shortcut_sheets)
        sheets = list(filter(is_background, sheets))
        for sheet in sheets:
            for pedigree in sheet.cohort.pedigrees:
                yield from pedigree.donors

    @listify
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample sheets.
        """
        for sub_step in self.sub_steps.values():
            if sub_step.name not in (LinkInStep.name,):
                yield from sub_step.get_result_files()

    def pick_kits_and_donors(self):
        """Return ``(library_kits, donors)`` with the donors with a matching kit and the kits with a
        matching donor.
        """
        kit_counts = {name: 0 for name in self.ngs_library_to_kit.values()}
        for name in self.ngs_library_to_kit.values():
            kit_counts[name] += 1
        donors = [
            donor
            for donor in self.all_donors()
            if donor.dna_ngs_library and donor.dna_ngs_library.name in self.ngs_library_to_kit
        ]
        return list(sorted(set(self.ngs_library_to_kit.values()))), donors, kit_counts

    def check_config(self):
        """Check that the necessary configuration is available for the step"""
        self.ensure_w_config(
            config_keys=("step_config", "sv_calling_targeted", "path_ngs_mapping"),
            msg="Path to NGS mapping not configured but required for targeted seq. CNV calling",
        )

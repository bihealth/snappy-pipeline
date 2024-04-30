"""SV calling for targeted sequencing

Based on the output of ``ngs_mapping``, call structural variants from depth of coverage,
read pair, and split read signal.
"""

import re

from biomedsheets.shortcuts import GermlineCaseSheet, is_not_background

from snappy_pipeline.utils import DictQuery, dictify, listify
from snappy_pipeline.workflows.abstract import BaseStep, WritePedigreeStepPart
from snappy_pipeline.workflows.common.delly import Delly2StepPart
from snappy_pipeline.workflows.common.gcnv.gcnv_run import RunGcnvStepPart
from snappy_pipeline.workflows.common.manta import MantaStepPart
from snappy_pipeline.workflows.common.melt import MeltStepPart
from snappy_pipeline.workflows.ngs_mapping import NgsMappingWorkflow

__author__ = "Manuel Holtgrewe <manuel.holtgrewe@bih-charite.de>"

#: Extensions of files to create as main payload (VCF)
EXT_VALUES = (".vcf.gz", ".vcf.gz.tbi", ".vcf.gz.md5", ".vcf.gz.tbi.md5")

#: Names of the files to create for the extension
EXT_NAMES = ("vcf", "vcf_tbi", "vcf_md5", "vcf_tbi_md5")

#: Available SV callers
SV_CALLERS = ("gcnv", "delly2", "manta", "melt")

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
      # Path to interval block list with PAR region for contig calling.
      path_par_intervals: null  # REQUIRED
      # Path to gCNV model - will execute analysis in CASE MODE.
      #
      # Example:
      #
      # - library: "Agilent SureSelect Human All Exon V6"  # Kit name, match in path_target_interval_list_mapping
      #   contig_ploidy: /path/to/ploidy-model         # Output from `DetermineGermlineContigPloidy`
      #   model_pattern: /path/to/model_*              # Output from `GermlineCNVCaller`
      precomputed_model_paths: []
      # Skip processing of the following libraries.  If the library is in
      # family/pedigree then all of the family/pedigree will be skipped.
      skip_libraries: []

    delly2:
      path_exclude_tsv: null  # optional
      map_qual: 1
      geno_qual: 5
      qual_tra: 20
      mad_cutoff: 9
      # Skip processing of the following libraries.  If the library is in
      # family/pedigree then all of the family/pedigree will be skipped.
      skip_libraries: []

    manta:
      num_threads: 16
      # Skip processing of the following libraries.  If the library is in
      # family/pedigree then all of the family/pedigree will be skipped.
      skip_libraries: []

    melt:
      me_refs_infix: 1KGP_Hg19
      me_types:
      - ALU
      - LINE1
      - SVA
      jar_file: REQUIRED
      genes_file: add_bed_files/1KGP_Hg19/hg19.genes.bed  # adjust, e.g., Hg38/Hg38.genes.bed
      # Skip processing of the following libraries.  If the library is in
      # family/pedigree then all of the family/pedigree will be skipped.
      skip_libraries: []
"""


class GcnvTargetedStepPart(RunGcnvStepPart):
    """Targeted seq. CNV calling with GATK4 gCNV"""

    def __init__(self, parent):
        super().__init__(parent)
        # Take shortcut from library to library kit.
        self.ngs_library_to_kit = self.parent.ngs_library_to_kit

    def get_params(self, action: str):
        param_fn = super().get_params(action)

        def _get_params(*args, **kwargs):
            params = param_fn(*args, **kwargs)
            if action == "contig_ploidy":
                par_intervals = (
                    self.config.get("helper_gcnv_model_targeted", {})
                    .get("gcnv", {})
                    .get("path_par_intervals", "")
                )
                params.update({"par_intervals": par_intervals})
            return params

        return _get_params


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
            (
                WritePedigreeStepPart,
                GcnvTargetedStepPart,
                Delly2StepPart,
                MantaStepPart,
                MeltStepPart,
            )
        )
        # Register sub workflows
        self.register_sub_workflow("ngs_mapping", self.config["path_ngs_mapping"])
        # Build dictionary with sample count per library kit
        _, _, self.library_kit_counts_dict = self.pick_kits_and_donors()

    @dictify
    def _build_ngs_library_to_kit(self):
        config = DictQuery(self.w_config).get("step_config/sv_calling_targeted/gcnv")
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
    def get_result_files(self):
        """Return list of result files for the NGS mapping workflow

        We will process all NGS libraries of all test samples in all sample sheets.
        """
        for sub_step in self.sub_steps.values():
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

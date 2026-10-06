import enum

from pydantic import model_validator

from snappy_pipeline.models import SnappyStepModel


class ToolNgsMapping(enum.StrEnum):
    BWA = "bwa"
    BWA_MEM2 = "bwa_mem2"


class ToolSomaticVariantCalling(enum.StrEnum):
    MUTECT2 = "mutect2"


class ToolGermlineVariantCalling(enum.StrEnum):
    GATK4_HC = "gatk4_hc"


class ToolVariantAnnotation(enum.StrEnum):
    VEP = "vep"


class InputVariantType(enum.StrEnum):
    CALLING = "calling"
    ANNOTATION = "annotation"
    FILTRATION = "filtration"


class RenameCombine(enum.StrEnum):
    TUMOR = "tumor"
    GERMLINE = "germline"


class CombineVariants(SnappyStepModel):
    tool_ngs_mapping: ToolNgsMapping = ToolNgsMapping.BWA

    somatic_variant_type: InputVariantType = InputVariantType.CALLING
    path_somatic_variant: str = "../somatic_variant_calling"
    tool_somatic_variant_calling: ToolSomaticVariantCalling = ToolSomaticVariantCalling.MUTECT2
    tool_somatic_variant_annotation: ToolVariantAnnotation | None = None
    is_somatic_variant_filtered: bool = False
    use_all_transcripts_for_somatic_annotations: bool = False

    germline_variant_type: InputVariantType = InputVariantType.CALLING
    path_germline_variant: str = "../germline_variant_calling"
    tool_germline_variant_calling: ToolGermlineVariantCalling = ToolGermlineVariantCalling.GATK4_HC
    tool_germline_variant_annotation: ToolVariantAnnotation | None = None
    is_germline_variant_filtered: bool = False
    use_all_transcripts_for_germline_annotations: bool = False

    rename_combined: RenameCombine | None = None

    @model_validator(mode="after")
    def ensure_variant_type_compatible_with_other_options(self):
        match self.germline_variant_type:
            case InputVariantType.CALLING:
                if self.is_germline_variant_filtered or self.tool_germline_variant_annotation:
                    raise ValueError(
                        "When the germline variant type is 'calling', "
                        "'is_germline_variant_filtered' must be set to false, "
                        "and 'tool_germline_variant_annotation' must be left unset"
                    )
            case InputVariantType.ANNOTATION:
                if not self.tool_germline_variant_annotation:
                    raise ValueError(
                        "When the germline variant type is 'annotation', the 'tool_germline_variant_annotation' must be set"
                    )
            case InputVariantType.FILTRATION:
                if not self.is_germline_variant_filtered:
                    raise ValueError(
                        "When the germline variant type is 'filtration', 'is_germline_variant_filtered' must be true"
                    )
        match self.somatic_variant_type:
            case InputVariantType.CALLING:
                if self.is_somatic_variant_filtered or self.tool_somatic_variant_annotation:
                    raise ValueError(
                        "When the somatic variant type is 'calling', "
                        "'is_somatic_variant_filtered' must be set to false, "
                        "and 'tool_somatic_variant_annotation' must be left unset"
                    )
            case InputVariantType.ANNOTATION:
                if not self.tool_somatic_variant_annotation:
                    raise ValueError(
                        "When the somatic variant type is 'annotation', the 'tool_somatic_variant_annotation' must be set"
                    )
            case InputVariantType.FILTRATION:
                if not self.is_somatic_variant_filtered:
                    raise ValueError(
                        "When the somatic variant type is 'filtration', 'is_somatic_variant_filtered' must be true"
                    )
        return self

    @model_validator(mode="after")
    def ensure_annotation_enabled_for_use_all_transcript(self):
        if self.use_all_transcripts_for_germline_annotations:
            if (
                self.germline_variant_type != InputVariantType.ANNOTATION
                or not self.tool_germline_variant_annotation
            ):
                raise ValueError(
                    "Using all germline transcripts can only be done if "
                    "the germline annotation tool is set and if "
                    "the germline variant type is 'annotation'"
                )
        if self.use_all_transcripts_for_somatic_annotations:
            if (
                self.somatic_variant_type != InputVariantType.ANNOTATION
                or not self.tool_somatic_variant_annotation
            ):
                raise ValueError(
                    "Using all somatic transcripts can only be done if "
                    "the somatic annotation tool is set and if "
                    "the somatic variant type is 'annotation'"
                )
        return self

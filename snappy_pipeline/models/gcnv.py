from typing import Annotated

from pydantic import ConfigDict, Field

from snappy_pipeline.models import SnappyModel


class TargetIntervalEntry(SnappyModel):
    """
    The following will match both the stock IDT library kit and the ones
    with spike-ins seen fromr Yale genomics.  The path above would be
    mapped to the name "default".
      - name: IDT_xGen_V1_0
        pattern: "xGen Exome Research Panel V1\\.0*"
        path: "path/to/targets.bed"
    """

    name: Annotated[str, Field(examples=["IDT_xGen_V1_0"])]

    pattern: Annotated[str, Field(examples=["xGen Exome Research Panel V1\\.0*"])]

    path: Annotated[str, Field(examples=["path/to/targets.bed"])]


class PrecomputedModelEntry(SnappyModel):
    model_config = ConfigDict(protected_namespaces=())

    library: Annotated[str, Field(examples=["Agilent SureSelect Human All Exon V6"])]
    """Kit name, match in path_target_interval_list_mapping"""

    contig_ploidy: Annotated[str, Field(examples=["/path/to/ploidy-model"])]
    """Output from `DetermineGermlineContigPloidy`"""

    model_pattern: Annotated[str, Field(examples=["/path/to/model_*"])]
    """Output from `GermlineCNVCaller`"""

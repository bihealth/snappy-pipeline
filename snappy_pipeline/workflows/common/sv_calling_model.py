"""Field types shared by the models of ``sv_calling_{targeted,wgs}``"""

from typing import Annotated

from snappy_pipeline.workflows.abstract.protocol import (
    DataSignature,
    DataType,
    ExpectedPathSchema,
    ExpectedRegions,
)

#: ``depends_on`` field type for the BED file of the genes of MELT
GeneRegions = Annotated[
    str, DataSignature(DataType.REGIONS, frozenset({"genes"})), ExpectedPathSchema(ExpectedRegions)
]

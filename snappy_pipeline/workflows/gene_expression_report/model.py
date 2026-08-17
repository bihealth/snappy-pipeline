from typing import Annotated

from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.gene_expression_quantification.model import ExpectedExpression


class GeneExpressionReportDependsOn(SnappyModel):
    gene_expression_quantification: Annotated[
        str,
        DataSignature(DataType.EXPRESSION, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedExpression),
    ] = "gene_expression_quantification"


class GeneExpressionReport(SnappyStepModel):
    depends_on: GeneExpressionReportDependsOn = Field(default_factory=GeneExpressionReportDependsOn)

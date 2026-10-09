from typing import Annotated

from snappy_pipeline.models import SnappyModel, SnappyStepModel
from snappy_pipeline.workflows.abstract.protocol import DataSignature, DataType, ExpectedPathSchema
from snappy_pipeline.workflows.gene_expression_quantification.model import ExpectedExpression


class GeneExpressionReportDependsOn(SnappyModel):
    expression: Annotated[
        str,
        DataSignature(DataType.EXPRESSION, frozenset({"rna"})),
        ExpectedPathSchema(ExpectedExpression),
    ]


class GeneExpressionReport(SnappyStepModel):
    depends_on: GeneExpressionReportDependsOn

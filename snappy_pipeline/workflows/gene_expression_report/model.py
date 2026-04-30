from pydantic import Field

from snappy_pipeline.models import SnappyModel, SnappyStepModel


class GeneExpressionReportDependsOn(SnappyModel):
    gene_expression_quantification: str = "gene_expression_quantification"


class GeneExpressionReport(SnappyStepModel):
    depends_on: GeneExpressionReportDependsOn = Field(default_factory=GeneExpressionReportDependsOn)

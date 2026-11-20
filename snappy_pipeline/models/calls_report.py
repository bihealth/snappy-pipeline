"""The ``calls_report`` section of the steps that write VCFs."""

from snappy_pipeline.models import SnappyModel


class CallsReport(SnappyModel):
    """The calls of each VCF of the task as a datavzrd site for the report, and the same table as a result file.
    The columns follow what the VCF carries: varlociraptor's event probabilities and observations, or the genotypes, allele frequencies and depths of other callers."""

    page_size: int = 25
    """Rows per page"""

    max_in_memory_rows: int = 50000
    """Rows up to which a view sorts and filters in the browser; datavzrd pages a larger view without sorting"""

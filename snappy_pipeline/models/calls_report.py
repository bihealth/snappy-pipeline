"""The ``calls_report`` section of the steps that write VCFs."""

import enum

from snappy_pipeline.models import SnappyModel


class Linkout(enum.StrEnum):
    gnomad = "gnomad"
    clinvar = "clinvar"
    dbsnp = "dbsnp"
    ensembl = "ensembl"
    genecards = "genecards"
    omim = "omim"
    varlociraptor_inspect = "varlociraptor_inspect"


class CallsReport(SnappyModel):
    """The calls of each VCF of the task as a datavzrd site for the report, and the same table as a result file.
    The columns follow what the VCF carries: varlociraptor's event probabilities and observations, or the genotypes, allele frequencies and depths of other callers; the gene, consequence, impact and HGVS of annotated variants; the gnomAD frequency and ClinVar classification that mehari adds."""

    page_size: int = 25
    """Rows per page"""

    max_in_memory_rows: int = 50000
    """Rows up to which a view sorts and filters in the browser; datavzrd pages a larger view without sorting"""

    linkouts: list[Linkout] = list(Linkout)
    """Links of a variant (``gnomad``, ``clinvar``, ``dbsnp``, and ``varlociraptor_inspect``, which opens varlociraptor's observations of the call without any names) and of its gene (``ensembl``, ``genecards``, ``omim``); empty for none.
    gnomAD, ClinVar and dbSNP by position and Ensembl need the genome build, which an annotation task in the chain tags its variants with; without one they are left out"""

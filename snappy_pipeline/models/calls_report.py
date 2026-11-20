"""The ``calls_report`` section of the steps that write VCFs."""

import enum

from snappy_pipeline.models import SnappyModel, ToolModel


class Linkout(enum.StrEnum):
    gnomad = "gnomad"
    clinvar = "clinvar"
    dbsnp = "dbsnp"
    ensembl = "ensembl"
    genecards = "genecards"
    omim = "omim"
    varlociraptor_inspect = "varlociraptor_inspect"


class IgvViews(ToolModel):
    """IGV views (igv-reports) of the top calls of each view in one page,
    ``igv/top.html`` in the site, which the ``igv`` column of the calls table links to;
    ``extra_arguments`` go to ``create_report``.
    The BAMs come from ``depends_on.alignments`` (the members of the entity and their matched normals),
    the genome from ``depends_on.reference``."""

    top: int = 20
    """Calls per view that get an IGV view, from the top of the view's ranking;
    a call in several views counts once"""

    views: list[str] = []
    """Views whose top calls get an IGV view:
    names such as ``present`` or ``ranked``, which must exist, or patterns such as ``denovo_*``,
    which may match none; empty takes every view"""

    flanking: int = 300
    """Width of the window around a call in base pairs, half on each side"""

    max_size_mb: int = 50
    """Size limit of the IGV page; the lowest-ranked calls are dropped until it fits"""


class CallsReport(SnappyModel):
    """The calls of each VCF of the task as a datavzrd site for the report,
    and the same table as a result file.
    The columns follow what the VCF carries:
    varlociraptor's event probabilities and observations, or the genotypes,
    allele frequencies and depths of other callers;
    the gene, consequence, impact and HGVS of annotated variants;
    the gnomAD frequency and ClinVar classification that mehari adds."""

    page_size: int = 25
    """Rows per page"""

    max_in_memory_rows: int = 50000
    """Rows up to which a view sorts and filters in the browser;
    datavzrd pages a larger view without sorting"""

    linkouts: list[Linkout] = list(Linkout)
    """Links of a variant (``gnomad``, ``clinvar``, ``dbsnp``, and ``varlociraptor_inspect``,
    which opens varlociraptor's observations of the call without any names)
    and of its gene (``ensembl``, ``genecards``, ``omim``); empty for none.
    gnomAD, ClinVar and dbSNP by position and Ensembl need the genome build,
    which an annotation task in the chain tags its variants with; without one they are left out"""

    igv: IgvViews | None = None
    """IGV views of the top calls of each view,
    from the BAMs of ``depends_on.alignments`` and the genome of ``depends_on.reference``; unset: none"""

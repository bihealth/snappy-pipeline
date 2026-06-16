"""Generic library-selection mixin for snappy pipeline step models.

Steps that need to enumerate a subset of NGS libraries from their samplesheet
inherit :class:`LibrarySelectionMixin`.  The selection is expressed as a
pandas ``DataFrame.query()`` string evaluated against the tidy library
DataFrame built by :func:`snappy_pipeline.workflows.abstract.build_library_dataframe`.
"""

from snappy_pipeline.models import SnappyModel


class LibrarySelectionMixin(SnappyModel):
    """Mixin for step models that support selecting a subset of libraries.

    Inherit this *before* ``SnappyStepModel`` so Pydantic sees the field:

    .. code-block:: python

        class MyStep(LibrarySelectionMixin, SnappyStepModel):
            ...
    """

    library_selection: str | None = None
    """Optional pandas ``DataFrame.query()`` expression to choose which
    libraries this task processes.

    The expression is evaluated against a tidy DataFrame with the following
    columns:

    ==================  ========================================================
    Column              Description
    ==================  ========================================================
    ``library_name``    NGS library identifier
    ``extraction_type`` ``"dna"``, ``"rna"``, … (always lower-cased)
    ``kind``            Study type: ``"cancer"`` or ``"germline"``
    ``role``            Biological/clinical function of the sample (see below)
    ``is_primary``      ``True`` for the primary tumor (cancer) or pedigree
                        index/proband (germline); ``False`` otherwise
    ``sex``             ``"male"``, ``"female"``, or ``"unknown"``
    ``donor_name``      Patient / family identifier
    ``sample_name``     Bio-sample name
    ==================  ========================================================

    **Role values**

    *Cancer sheets:*
      ``"tumor"`` — ``isTumor=True`` samples |br|
      ``"normal"`` — ``isTumor=False`` samples

    *Germline / trio sheets:*
      ``"index"`` — index / proband of the pedigree |br|
      ``"father"`` — paternal donor (has ``father_pk`` defined) |br|
      ``"mother"`` — maternal donor (has ``mother_pk`` defined) |br|
      ``"affected"`` — affected non-index members |br|
      ``"unaffected"`` — unaffected members

    **Examples**

    .. code-block:: yaml

        # Only primary tumor DNA libraries (default for cancer sheets)
        library_selection: "role == 'tumor' and extraction_type == 'dna'"

        # Only normal (paired-normal) libraries
        library_selection: "role == 'normal' and extraction_type == 'dna'"

        # Primary tumor only (exclude relapse/additional tumors)
        library_selection: "role == 'tumor' and is_primary"

        # Index proband only (germline/trio)
        library_selection: "role == 'index' and extraction_type == 'dna'"

        # Tumor OR relapse samples
        library_selection: "role in ('tumor', 'relapse')"

        # Explicit allowlist
        library_selection: "library_name in ['lib-001', 'lib-002']"

    When ``None``, a per-sheet-type default is applied automatically:

    * ``"cancer"``   →  ``"role == 'tumor' and extraction_type == 'dna'"``
    * ``"germline"`` →  ``"extraction_type == 'dna'"``
    """  # noqa: E501

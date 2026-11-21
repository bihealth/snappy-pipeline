"""The strand protocol of RNA libraries.

RSeQC's ``infer_experiment`` decides it from the reads on the exons of ``depends_on.exons``, unless the configuration forces it.
``ngs_mapping`` with ``tool: star`` and ``gene_expression_quantification`` with ``tool: strandedness`` run it through the ``rseqc`` wrapper, which writes the decision as JSON.
"""

import enum
from typing import Any

from snappy_pipeline.models import ToolModel


class Strand(enum.IntEnum):
    """The strand protocol of a library, numbered like the counts columns of STAR and the ``-s`` of featureCounts"""

    infer = -1
    unstranded = 0
    forward = 1
    reverse = 2


class Strandedness(ToolModel):
    """The strand decision of ``infer_experiment``; ``extra_arguments`` go to ``infer_experiment.py``"""

    strand: Strand = Strand.infer
    """-1 infers the strand, 0, 1 and 2 force unstranded, forward and reverse"""

    threshold: float = 0.85
    """Minimum proportion of reads mapped to forward/reverse direction to call the protocol"""


def strandedness_params(section: Strandedness, library_name: str) -> dict[str, Any]:
    """Return the params of the ``rseqc`` wrapper for one library."""
    return {"config": section.model_dump(by_alias=True), "library_name": library_name}

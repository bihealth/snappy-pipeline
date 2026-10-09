from collections.abc import Mapping

from snappy_pipeline.models import SnappyStepModel


def require_tool_dependencies(model: SnappyStepModel, required: Mapping[str, tuple[str, ...]]):
    """Raise ``ValueError`` unless every ``depends_on`` field that ``model.tool`` reads is set.

    ``required`` maps a tool to the ``depends_on`` fields it reads.
    """
    missing = [
        field for field in required.get(model.tool, ()) if not getattr(model.depends_on, field)
    ]
    if missing:
        fields = " and ".join(f"depends_on.{field}" for field in missing)
        raise ValueError(f"tool={model.tool} needs {fields}")

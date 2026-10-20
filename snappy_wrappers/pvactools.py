"""Helpers shared by the pVACtools wrappers."""

import os


def bindings_for_container(
    files_to_bind,
    base_container_path: str,
    access: str = "ro",
    exclude_bind: list[str] | tuple[str, ...] = (),
) -> tuple[str, dict[str, str]]:
    """Bind the directories of the (named, non-list) files into a container.

    Paths are resolved with ``realpath``, and each directory is bound once, as
    ``<base_container_path>/d<N>``.
    Returns the ``--bind`` options and the path of each file inside the container.
    Keys in ``exclude_bind`` are skipped.
    """
    fns = {}
    for key, path in files_to_bind.items():
        if key in exclude_bind:
            continue
        assert not isinstance(path, list), f"Input files cannot be lists (offender is '{key}')"
        assert key not in fns, f"Duplicate input '{key}'"
        fns[key] = os.path.realpath(path)

    dirs = {key: os.path.dirname(path) for key, path in fns.items()}
    bindings = {d: i for i, d in enumerate(set(dirs.values()))}

    container_fns = {
        key: f"{base_container_path}/d{bindings[dirs[key]]}/{os.path.basename(path)}"
        for key, path in fns.items()
    }
    binding_cmd = " ".join(
        f"--bind {d}:{base_container_path}/d{i}:{access}" for d, i in bindings.items()
    )
    return (binding_cmd, container_fns)

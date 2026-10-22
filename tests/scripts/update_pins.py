#!/usr/bin/env python3
"""Pin wrapper environments again with snakedeploy (needs conda).

Without arguments, every environment file that has a pin file; a symlinked environment file
shares its target's pin. With arguments, the given environment files, also new ones.

Usage: update_pins.py [environment.yaml ...]
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

WRAPPERS = Path(__file__).resolve().parents[2] / "snappy_wrappers" / "wrappers"
PIN_FILE = "environment.linux-64.pin.txt"


def pinned_environments() -> list[Path]:
    """Return the environment files that have their own pin file."""
    return sorted(
        env
        for env in WRAPPERS.rglob("environment.yaml")
        if not env.is_symlink() and (env.parent / PIN_FILE).is_file()
    )


def main() -> int:
    envs = sys.argv[1:] or [str(env) for env in pinned_environments()]
    cmd = ["snakedeploy", "pin-conda-envs", "--conda-frontend", "conda", *envs]
    return subprocess.run(cmd, check=False).returncode


if __name__ == "__main__":
    sys.exit(main())

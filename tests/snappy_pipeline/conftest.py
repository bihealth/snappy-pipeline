"""Fixtures shared by the tests of ``snappy_pipeline``."""

import os
import subprocess

import pytest

#: A fixed commit identity, and no configuration of the user or the system
GIT_ENV = {
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_NOSYSTEM": "1",
    "GIT_AUTHOR_NAME": "snappy",
    "GIT_AUTHOR_EMAIL": "snappy@example.org",
    "GIT_COMMITTER_NAME": "snappy",
    "GIT_COMMITTER_EMAIL": "snappy@example.org",
}


@pytest.fixture
def git(monkeypatch):
    """Return ``git(directory, *args)``, which runs git there with :data:`GIT_ENV` and returns
    its output; the environment applies to the test as well."""
    for key, value in GIT_ENV.items():
        monkeypatch.setenv(key, value)

    def run(directory, *args) -> str:
        return subprocess.run(
            ["git", "-C", str(directory), *args], capture_output=True, text=True, check=True
        ).stdout

    return run

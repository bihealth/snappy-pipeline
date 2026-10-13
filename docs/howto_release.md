(howto_release)=

# How To: Release

Releases are made by release-please (`.github/workflows/release-please.yml`).
It reads the Conventional Commit messages on `main`, so every squash commit needs a correct type (`feat`, `fix`, ...).
Breaking changes use `!` or a `BREAKING CHANGE` footer.

1. Merge pull requests into `main`. release-please opens or updates a release pull request that bumps `snappy_pipeline/version.py` and adds the entries to `CHANGELOG.md`.
2. Review the release pull request and merge it. release-please then creates the tag `v<version>` and the GitHub release.

Packages are not built and the documentation is not deployed by this workflow.

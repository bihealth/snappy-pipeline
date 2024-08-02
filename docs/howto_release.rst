.. _howto_release:

===============
How To: Release
===============

1. Update the version in the following files:

    - ``installation.rst`` (look for ``snappy_pipeline.git@VERSION``)
    - **TODO** more?

2. Create a tag and push it

    .. code-block:: shell

        $ git tag v0.1.0
        $ git push --tags origin

That's it, so far we don't create packages or deploy the documentation.

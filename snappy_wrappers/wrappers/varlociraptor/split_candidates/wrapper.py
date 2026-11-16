# -*- coding: utf-8 -*-
"""Split the candidates of an entity into chunks of about equal work for varlociraptor.

The work of ``preprocess`` grows with the reads at each candidate, so a candidate weighs its
``INFO/DP`` (one without it weighs one read), and every chunk gets about the same summed depth.
The records of a breakend event (``INFO/EVENT``, or mates by ``INFO/MATEID``) stay in one chunk,
as varlociraptor evaluates them together. The chunks keep the input order, so they are sorted.
"""

import pysam


def info(record, key: str, defined: set[str]):
    """The INFO value of a record, ``None`` when the header does not define the key (pysam
    raises for undefined keys) or the record lacks it."""
    return record.info.get(key) if key in defined else None


def weight(record, defined: set[str]) -> int:
    """The depth of a candidate, one for a record without ``INFO/DP``."""
    depth = info(record, "DP", defined)
    if isinstance(depth, tuple):
        depth = depth[0] if depth else None
    return max(int(depth), 1) if depth else 1


def group(record, defined: set[str]):
    """The key of the breakend event of a record, ``None`` for a record evaluated alone."""
    event = info(record, "EVENT", defined)
    if event:
        return ("event", event if isinstance(event, str) else event[0])
    mates = info(record, "MATEID", defined)
    if mates:
        mates = (mates,) if isinstance(mates, str) else tuple(mates)
        return ("mates", tuple(sorted({*mates, record.id or ""})))
    return None


def plan(weights, groups, num_chunks: int) -> list[int]:
    """Return the chunk of each record: by cumulative weight, the records of a group in the
    chunk of its first record."""
    total = sum(weights)
    quota = total / num_chunks if total else 1
    chunk_of_group = {}
    chunks = []
    done = 0
    for record_weight, key in zip(weights, groups):
        chunk = min(int(done / quota), num_chunks - 1)
        if key is not None:
            chunk = chunk_of_group.setdefault(key, chunk)
        chunks.append(chunk)
        done += record_weight
    return chunks


def main(candidates: str, outputs: list[str]) -> None:
    with pysam.VariantFile(candidates) as reader:
        defined = set(reader.header.info.keys())
        weights, groups = [], []
        for record in reader:
            weights.append(weight(record, defined))
            groups.append(group(record, defined))
    chunks = plan(weights, groups, len(outputs))
    with pysam.VariantFile(candidates) as reader:
        writers = [pysam.VariantFile(path, "wb", header=reader.header) for path in outputs]
        try:
            for record, chunk in zip(reader, chunks):
                writers[chunk].write(record)
        finally:
            for writer in writers:
                writer.close()


if __name__ == "__main__":
    if snakemake := locals().get("snakemake", None):
        from snappy_wrappers.snappy_wrapper import PythonWrapper

        PythonWrapper(snakemake).run(
            lambda: main(snakemake.input.candidates, list(snakemake.output.candidates))
        )

"""Varlociraptor scenarios rendered from the library table.

A scenario names the samples of an entity (their priors or allele frequency universes) and the
events to compute probabilities for. The sample keys are library names, so the VCF samples are
named like everywhere in snappy. The presets ``germline`` and ``somatic`` are built here; a
user's YTE template replaces the preset's YAML.
"""

import re
from typing import Any

import yaml

#: What varlociraptor accepts as a sample or event name
IDENTIFIER = re.compile(r"^[A-Za-z0-9_.-]+$")

#: The sample key of the unobserved normal of a tumor without a matched normal
VIRTUAL_NORMAL = "normal"

#: The library table's values for a parent that is unknown or has no DNA library
NO_PARENT = ("0", "", None)


def is_ffpe(value: Any) -> bool:
    """Whether a sample sheet value of ``isFfpe`` (``Y``/``N``, ``true``/``false``) marks an FFPE
    library; biomedsheets keeps custom field strings, so the words are read here."""
    return str(value).strip().lower() in ("y", "yes", "true", "1")


def germline_scenario(samples: list[dict[str, Any]], species: dict[str, Any]) -> dict[str, Any]:
    """Return the ``germline`` preset for the library rows ``samples`` of a pedigree.

    Every member is a sample with its ``sex``, or diploid everywhere when the sex is unknown; a
    member whose parents are both members inherits from them (mendelian). The events are
    ``present`` (a member carries the variant) and ``denovo`` (a child carries a variant that
    its parents do not), with ``present`` excluding it. Without a member with both parents,
    ``denovo`` is the impossible event ``false``, so one FDR control fits every pedigree.

    Varlociraptor models inheritance from both parents only, so a member with a parent outside
    the samples (one parent sequenced, or a parent that the task does not take) raises
    ``ValueError``: the preset would otherwise treat related libraries as unrelated.
    """
    names = [sample["library_name"] for sample in samples]
    scenario: dict[str, Any] = {
        "species": {
            "heterozygosity": species["heterozygosity"],
            "germline-mutation-rate": species["germline_mutation_rate"],
            "ploidy": species["ploidy"],
        },
        "samples": {},
        "events": {},
    }
    children = []
    for sample in samples:
        name = sample["library_name"]
        entry: dict[str, Any] = {}
        if sample.get("sex") in ("male", "female"):
            entry["sex"] = sample["sex"]
        else:
            entry["ploidy"] = 2
        if is_ffpe(sample.get("is_ffpe")):
            entry["fragment_end_bias"] = True
        father, mother = sample.get("father_name"), sample.get("mother_name")
        known = [parent for parent in (father, mother) if parent not in NO_PARENT]
        present = [parent for parent in known if parent in names]
        if len(present) == 2:
            entry["inheritance"] = {"mendelian": {"from": [mother, father]}}
            children.append((name, mother, father))
        elif known:
            raise ValueError(
                f"{name} has the parents {' and '.join(known)} in the sample sheet, but the "
                "germline scenario models inheritance only from both parents in the task "
                f"({', '.join(present) or 'none of them'} here); take both parents into the "
                "task or give a scenario template"
            )
        scenario["samples"][name] = entry
    if children:
        # A formula is a conjunction, a disjunction or an atom; only the terms of a disjunction
        # may be parenthesised conjunctions
        terms = [
            f"{child}:]0.0,1.0] & {mother}:0.0 & {father}:0.0" for child, mother, father in children
        ]
        scenario["events"]["denovo"] = (
            terms[0] if len(terms) == 1 else " | ".join(f"({term})" for term in terms)
        )
        scenario["events"]["present"] = "!($denovo | $absent)"
    else:
        # The literal false may not be negated, so present names the members
        scenario["events"]["denovo"] = False
        scenario["events"]["present"] = " | ".join(f"{name}:]0.0,1.0]" for name in names)
    return scenario


def somatic_scenario(
    tumor: str, normal: str | None, purity: float | None, ffpe: frozenset[str] = frozenset()
) -> dict[str, Any]:
    """Return the ``somatic`` preset of dna-seq-varlociraptor for a tumor and its matched normal.

    Without a matched normal, the normal stays in the scenario as an unobserved sample, so the
    events are the same for every tumor. A known purity is the tumor's contamination by the
    normal; without one, the tumor counts as pure. The libraries in ``ffpe`` get varlociraptor's
    model of the artifacts of an FFPE fragment (``fragment_end_bias``).
    """
    normal = normal or VIRTUAL_NORMAL
    tumor_sample: dict[str, Any] = {"resolution": 0.01, "universe": "[0.0,1.0]"}
    if purity is not None:
        tumor_sample["contamination"] = {"by": normal, "fraction": round(1.0 - purity, 6)}
    normal_sample: dict[str, Any] = {"resolution": 0.1, "universe": "0.0 | 0.5 | 1.0 | ]0.0,0.5["}
    for name, sample in ((tumor, tumor_sample), (normal, normal_sample)):
        if name in ffpe:
            sample["fragment_end_bias"] = True
    return {
        "samples": {tumor: tumor_sample, normal: normal_sample},
        "events": {
            "somatic_tumor_low": f"{tumor}:]0.0,0.1[ & {normal}:0.0",
            "somatic_tumor_medium": f"{tumor}:[0.1,0.3[ & {normal}:0.0",
            "somatic_tumor_high": f"{tumor}:[0.3,1.0] & {normal}:0.0",
            "somatic_normal": f"{normal}:]0.0,0.5[",
            "germline_hom": f"{normal}:1.0",
            "germline_het": f"{normal}:0.5",
        },
    }


def render_template(
    template: str, samples: list[dict[str, Any]], library_name: str, settings: dict[str, Any]
) -> dict[str, Any]:
    """Return the scenario of a YTE template, rendered with the variables ``samples`` (the
    library rows of the entity's libraries), ``library_name`` (the entity) and ``settings``."""
    import yte

    variables = {"samples": samples, "library_name": library_name, "settings": settings}
    scenario = yte.process_yaml(template, variables=variables)
    if not isinstance(scenario, dict) or not {"samples", "events"} <= set(scenario):
        raise ValueError(f"The scenario template of {library_name} needs samples and events")
    return scenario


def check_scenario(scenario: dict[str, Any], members: list[str]) -> None:
    """Raise ``ValueError`` unless every observed library is a sample of the scenario."""
    missing = sorted(set(members) - set(scenario["samples"]))
    if missing:
        raise ValueError(
            f"The scenario lacks the samples {', '.join(missing)}; its samples: "
            + ", ".join(scenario["samples"])
        )


def to_yaml(scenario: dict[str, Any]) -> str:
    """Return the scenario as YAML for ``varlociraptor call variants generic --scenario``."""
    return yaml.safe_dump(scenario, sort_keys=False, width=float("inf"))

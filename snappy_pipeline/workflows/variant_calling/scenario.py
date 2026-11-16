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


def germline_scenario(samples: list[dict[str, Any]], species: dict[str, Any]) -> dict[str, Any]:
    """Return the ``germline`` preset for the library rows ``samples`` of a pedigree.

    Every member is a sample with its ``sex``, or diploid everywhere when the sex is unknown; a
    member whose parents are both members inherits from them (mendelian). The events are
    ``present`` (a member carries the variant) and, when a member has both parents, ``denovo``
    (a child carries a variant that its parents do not), with ``present`` then excluding it.
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
        father, mother = sample.get("father_name"), sample.get("mother_name")
        if father in names and mother in names:
            entry["inheritance"] = {"mendelian": {"from": [mother, father]}}
            children.append((name, mother, father))
        scenario["samples"][name] = entry
    if children:
        scenario["events"]["denovo"] = " | ".join(
            f"({child}:]0.0,1.0] & {mother}:0.0 & {father}:0.0)"
            for child, mother, father in children
        )
        scenario["events"]["present"] = "!($denovo | $absent)"
    else:
        scenario["events"]["present"] = " | ".join(f"{name}:]0.0,1.0]" for name in names)
    return scenario


def somatic_scenario(tumor: str, normal: str | None, purity: float | None) -> dict[str, Any]:
    """Return the ``somatic`` preset of dna-seq-varlociraptor for a tumor and its matched normal.

    Without a matched normal, the normal stays in the scenario as an unobserved sample, so the
    events are the same for every tumor. A known purity is the tumor's contamination by the
    normal; without one, the tumor counts as pure.
    """
    normal = normal or VIRTUAL_NORMAL
    tumor_sample: dict[str, Any] = {"resolution": 0.01, "universe": "[0.0,1.0]"}
    if purity is not None:
        tumor_sample["contamination"] = {"by": normal, "fraction": round(1.0 - purity, 6)}
    return {
        "samples": {
            tumor: tumor_sample,
            normal: {"resolution": 0.1, "universe": "0.0 | 0.5 | 1.0 | ]0.0,0.5["},
        },
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

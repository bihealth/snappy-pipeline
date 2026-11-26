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

#: The library table's ``disease_state`` (the PED phenotype) of affected and unaffected members
AFFECTED, UNAFFECTED = 2, 1


def is_ffpe(value: Any) -> bool:
    """Whether a sample sheet value of ``is_ffpe`` (``true``/``false``, or the older ``Y``/``N``) marks an FFPE
    library; the sheet column is a boolean, the words are read for older tables."""
    return str(value).strip().lower() in ("y", "yes", "true", "1")


def _excluding(formula: str, exclusions: list[str]) -> str:
    """Return ``formula`` with the negated ``exclusions`` (``$`` references) conjoined."""
    return " & ".join([formula, *(f"!{excluded}" for excluded in exclusions)])


def event_suffix(sample: dict[str, Any]) -> str:
    """The member's part of a per-member event name such as ``denovo_<member>``:
    the donor name (else the library name) in lower case with letters, digits and ``_`` only,
    so that ``PROB_<EVENT>`` is a valid VCF INFO ID and the name comes back from the header unchanged."""
    name = str(sample.get("donor_name") or sample["library_name"]).lower()
    return re.sub(r"[^a-z0-9]", "_", name)


def germline_scenario(samples: list[dict[str, Any]], species: dict[str, Any]) -> dict[str, Any]:
    """Return the ``germline`` preset for the library rows ``samples`` of a pedigree.

    Every sample has its ``sex``, or is diploid everywhere when the sex is unknown.
    A sample whose parents are both samples inherits from them (mendelian).
    A sample with ``observed: False`` is an unsequenced relative that connects members,
    or the unknown parent of a member;
    its ``library_name`` is its donor name.
    It has no observations, so the posteriors integrate over its genotype,
    and inheritance flows through it.

    The events partition the configurations in which an observed member carries the variant,
    by the sheet's ``isAffected``.
    ``disease_state`` 2 is affected, 1 unaffected, anything else unknown, which constrains nothing.
    ``denovo_<member>``, one per observed child whose parents are both samples:
    the child carries a variant that its parents do not (a de novo mutation in that member,
    which its descendants inherit).
    ``<member>`` is the donor name in lower case (:func:`event_suffix`).
    ``recessive``: every affected member is homozygous, no unaffected member is (``!member:1.0``).
    ``dominant``: every affected member carries the variant, no unaffected member does.
    ``present``: the rest.
    Each event excludes the ones before it, so they are disjoint.
    ``dominant`` thus means inherited dominant; in a trio with unaffected parents it is empty.
    Without an affected member, ``recessive`` and ``dominant`` are the impossible event ``false``,
    and a pedigree without a child has no de novo event,
    so FDR control over ``denovo_*``, ``recessive``,
    ``dominant`` and ``present`` fits every pedigree.

    The patterns of the events are ``expressions``,
    so that the events negate ``$`` references without a negation inside.
    varlociraptor 8.12.0 expands references one level deep and panics on a negated reference to an event that negates another.

    Varlociraptor models inheritance from both parents only.
    A member with a parent outside the samples (one parent sequenced,
    or a parent the task does not take) raises ``ValueError``;
    the preset would otherwise treat related libraries as unrelated.
    The step adds unsequenced parents as unobserved samples first (``unknown_parents``).
    """
    names = [sample["library_name"] for sample in samples]
    observed = [sample["library_name"] for sample in samples if sample.get("observed", True)]
    affected = [name for name in observed if _disease_state(samples, name) == AFFECTED]
    unaffected = [name for name in observed if _disease_state(samples, name) == UNAFFECTED]
    scenario: dict[str, Any] = {
        "species": {
            "heterozygosity": species["heterozygosity"],
            "germline-mutation-rate": species["germline_mutation_rate"],
            "ploidy": species["ploidy"],
        },
        "samples": {},
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
            if name in observed:
                children.append((name, mother, father, event_suffix(sample)))
        elif known:
            raise ValueError(
                f"{name} has the parents {' and '.join(known)} in the sample sheet, but the "
                "germline scenario models inheritance only from both parents in the task "
                f"({', '.join(present) or 'none of them'} here); take both parents into the "
                "task or give a scenario template"
            )
        scenario["samples"][name] = entry
    expressions: dict[str, str] = {}
    events: dict[str, Any] = {}
    exclusions: list[str] = []
    named: dict[str, str] = {}
    for child, mother, father, suffix in children:
        event = f"denovo_{suffix}"
        if event in named:
            raise ValueError(
                f"{child} and {named[event]} both get the event {event}; give their donors "
                "names that differ in letters and digits"
            )
        named[event] = child
        pattern = f"{event}_pattern"
        expressions[pattern] = f"{child}:]0.0,1.0] & {mother}:0.0 & {father}:0.0"
        events[event] = _excluding(f"${pattern}", exclusions)
        exclusions.append(f"${pattern}")
    if affected:
        # "not homozygous" as a negated atom: varlociraptor 8.12.1 mishandles an interval that
        # holds 0.0 and 0.5, such as [0.0,1.0[, and the negation holds on every ploidy
        expressions["recessive_pattern"] = " & ".join(
            [f"{name}:1.0" for name in affected] + [f"!{name}:1.0" for name in unaffected]
        )
        events["recessive"] = _excluding("$recessive_pattern", exclusions)
        exclusions.append("$recessive_pattern")
        expressions["dominant_pattern"] = " & ".join(
            [f"{name}:]0.0,1.0]" for name in affected] + [f"{name}:0.0" for name in unaffected]
        )
        events["dominant"] = _excluding("$dominant_pattern", exclusions)
        exclusions.append("$dominant_pattern")
    else:
        events["recessive"] = events["dominant"] = False
    carried = " | ".join(f"{name}:]0.0,1.0]" for name in observed)
    events["present"] = _excluding(f"({carried})", exclusions) if exclusions else carried
    if expressions:
        scenario["expressions"] = expressions
    scenario["events"] = events
    return scenario


def _disease_state(samples: list[dict[str, Any]], name: str) -> Any:
    for sample in samples:
        if sample["library_name"] == name:
            return sample.get("disease_state")
    return None


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
    """Raise ``ValueError`` unless every observed library is a sample of the scenario and every
    sample name is an identifier that varlociraptor accepts."""
    missing = sorted(set(members) - set(scenario["samples"]))
    if missing:
        raise ValueError(
            f"The scenario lacks the samples {', '.join(missing)}; its samples: "
            + ", ".join(scenario["samples"])
        )
    invalid = sorted(name for name in scenario["samples"] if not IDENTIFIER.match(name))
    if invalid:
        raise ValueError(
            "varlociraptor needs sample names of letters, digits, '_', '.' and '-', unlike "
            + ", ".join(invalid)
        )


def to_yaml(scenario: dict[str, Any]) -> str:
    """Return the scenario as YAML for ``varlociraptor call variants generic --scenario``."""
    return yaml.safe_dump(scenario, sort_keys=False, width=float("inf"))

# -*- coding: utf-8 -*-
"""Pedigree manipulation related methods"""
from biomedsheets.shortcuts import KEY_IS_AFFECTED, KEY_SEX


def append_pedigree_to_ped(pedigree, output_path):
    """Append pedigree to ped file

    Based on method from ``biomedsheets``, adjusted to just out sample id instead of library name
    (e.g., 'P001' instead of 'P001-N1-DNA1-WGS1').

    :param pedigree: Pedigree.
    :type pedigree: biomedsheets.shortcuts.germline.Pedigree

    :param output_path: Path to output ped file.
    :type output_path: str
    """
    with open(output_path, "wt") as f:
        family = "FAM_" + pedigree.index.name
        for donor in pedigree.donors:
            # Parse affected status and sex
            affected = {"affected": "2", "unaffected": "1", "unknown": "0"}[
                donor.extra_infos.get(KEY_IS_AFFECTED, "unknown")
            ]
            sex = {"male": "1", "female": "2", "unknown": "0"}[
                donor.extra_infos.get(KEY_SEX, "unknown")
            ]

            # Parse father sample name
            father = "0"
            if donor.father_pk:
                if hasattr(pedigree.pk_to_donor[donor.father_pk], "dna_ngs_library"):
                    donor_father = pedigree.pk_to_donor[donor.father_pk]
                    father = donor_father.name

            # Parse mother sample name
            mother = "0"
            if donor.mother_pk:
                if hasattr(pedigree.pk_to_donor[donor.mother_pk], "dna_ngs_library"):
                    donor_mother = pedigree.pk_to_donor[donor.mother_pk]
                    mother = donor_mother.name

            # Parse sample name and output to file
            if hasattr(donor, "dna_ngs_library"):
                name = donor.name
                print("\t".join((family, name, father, mother, sex, affected)), file=f)

"""
Donor quality analysis for coordination shells.

Detects bidentate chelation pairs (same residue contributing 2+ donors)
and classifies individual donor atom types. Bidentate chelation count
is the strongest single predictor of lanthanide-binding peptide affinity
(r = −0.67 vs log K_D), motivated by the thermodynamic chelate effect.

Reference: Gutenthaler-Tietze et al. (2025) — retro EF-hand peptides.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from terrahedral.core import CoordinationSite


# ── Donor type classification ──────────────────────────────────────────────

# Canonical donor types and their display colours (CSS variable fallbacks)
DONOR_COLORS = {
    "carboxylate":  "#e05050",  # red — Asp/Glu COO⁻
    "backbone":     "#e0a050",  # orange — backbone C=O
    "water":        "#70a0e0",  # blue — coordinated water
    "amide":        "#a070e0",  # purple — Asn/Gln amide
    "hydroxyl":     "#50c070",  # green — Ser/Thr –OH
    "thiolate":     "#e0e050",  # yellow — Cys S⁻
    "imidazole":    "#50b0b0",  # teal — His
    "phenolate":    "#e070a0",  # pink — Tyr
    "phosphate":    "#c08040",  # brown — PO₄
    "unknown":      "#909090",  # grey
}

DONOR_LABELS = {
    "carboxylate":  "COO⁻",
    "backbone":     "C=O",
    "water":        "H₂O",
    "amide":        "amide",
    "hydroxyl":     "–OH",
    "thiolate":     "S⁻",
    "imidazole":    "His",
    "phenolate":    "Tyr",
    "phosphate":    "PO₄",
    "unknown":      "other",
}


def classify_donor(residue: str, atom_name: str, element: str = "O") -> str:
    """
    Classify a ligand donor atom by its chemical type.

    Parameters
    ----------
    residue : residue name (e.g. "ASP", "HOH", "WAT")
    atom_name : atom name within residue (e.g. "OD1", "O", "OE2")
    element : element symbol

    Returns
    -------
    One of the keys in DONOR_COLORS.
    """
    res = (residue or "").upper().strip()
    aname = (atom_name or "").upper().strip()

    # Water
    if any(tag in res for tag in ("WAT", "HOH", "H2O", "DOD")):
        return "water"

    # Phosphate
    if "PO4" in res or "HPO" in res:
        return "phosphate"

    # Protein sidechain donors
    if ("ASP" in res and aname in ("OD1", "OD2")):
        return "carboxylate"
    if ("GLU" in res and aname in ("OE1", "OE2")):
        return "carboxylate"
    if ("ASN" in res and aname in ("OD1", "ND2")):
        return "amide"
    if ("GLN" in res and aname in ("OE1", "NE2")):
        return "amide"
    if ("CYS" in res and aname == "SG"):
        return "thiolate"
    if ("HIS" in res and aname in ("ND1", "NE2")):
        return "imidazole"
    if ("TYR" in res and aname == "OH"):
        return "phenolate"
    if (("THR" in res or "SER" in res) and aname in ("OG", "OG1")):
        return "hydroxyl"

    # Backbone carbonyl
    if aname == "O":
        return "backbone"

    return "unknown"


# ── Bidentate chelation detection ──────────────────────────────────────────

CHELATION_COLORS = {
    "succinate-like":         "#9060d0",   # purple — backbone + sidechain
    "bidentate-carboxylate":  "#50b0b0",   # teal — both sidechain COO⁻
    "bidentate-other":        "#909090",   # grey
}


def find_bidentate_pairs(site: "CoordinationSite") -> list[dict]:
    """
    Find residues contributing 2+ donor atoms to the metal centre.

    Each such pair creates a chelate ring, which thermodynamically
    favours binding via the chelate effect (entropic stabilisation).

    Returns
    -------
    List of dicts, each with keys:
        residue     : residue identifier string
        donors      : list of {label, element, atom_name, distance, donor_type}
        type        : 'succinate-like' | 'bidentate-carboxylate' | 'bidentate-other'
        mean_distance : mean M–L distance for donors in this pair
    """
    bd = site.bond_dict

    # Group donors by residue
    residue_donors: dict[str, list[dict]] = {}
    for lig in site.ligands:
        res = lig.residue or ""
        # Skip waters — they can't form bidentate chelates
        res_upper = res.upper()
        if any(tag in res_upper for tag in ("WAT", "HOH", "H2O", "DOD")):
            continue
        if not res:
            continue

        dist = bd.get(lig.label, lig.distance_to(site.metal))
        dtype = classify_donor(res, lig.atom_name or lig.label, lig.element)

        if res not in residue_donors:
            residue_donors[res] = []
        residue_donors[res].append({
            "label": lig.label,
            "element": lig.element,
            "atom_name": lig.atom_name or "",
            "distance": dist,
            "donor_type": dtype,
        })

    pairs = []
    for res, donors in residue_donors.items():
        if len(donors) < 2:
            continue

        atom_names = [d["atom_name"] for d in donors]
        has_backbone = "O" in atom_names  # backbone carbonyl
        has_sidechain = any(
            a in ("OD1", "OD2", "OE1", "OE2") for a in atom_names
        )

        if has_backbone and has_sidechain:
            chel_type = "succinate-like"
        elif all(a in ("OD1", "OD2", "OE1", "OE2") for a in atom_names):
            chel_type = "bidentate-carboxylate"
        else:
            chel_type = "bidentate-other"

        mean_d = sum(d["distance"] for d in donors) / len(donors)
        pairs.append({
            "residue": res,
            "donors": donors,
            "type": chel_type,
            "mean_distance": round(mean_d, 3),
        })

    return pairs


# ── Full donor analysis for a site ─────────────────────────────────────────

def analyze_donors(site: "CoordinationSite") -> dict:
    """
    Complete donor analysis for a coordination site.

    Returns dict with:
        donor_types   : list of type strings, one per ligand
        donor_summary : {type: count} for display
        n_protein     : number of protein-derived donors (non-water)
        n_water       : number of water donors
        bidentate     : list from find_bidentate_pairs()
        n_bidentate   : count of bidentate pairs
    """
    bd = site.bond_dict
    donor_types = []
    summary: dict[str, int] = {}

    for lig in site.ligands:
        dtype = classify_donor(
            lig.residue or "",
            lig.atom_name or lig.label,
            lig.element,
        )
        donor_types.append(dtype)
        summary[dtype] = summary.get(dtype, 0) + 1

    n_water = summary.get("water", 0)
    n_protein = len(donor_types) - n_water

    bidentate = find_bidentate_pairs(site)

    return {
        "donor_types": donor_types,
        "donor_summary": summary,
        "n_protein": n_protein,
        "n_water": n_water,
        "bidentate": bidentate,
        "n_bidentate": len(bidentate),
    }

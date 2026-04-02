"""
Geometry classification for coordination sites.

Implements Addison τ₅ (5-coordinate), Yang τ₄ (4-coordinate),
and heuristic geometry assignment.
"""

from __future__ import annotations

from typing import Optional

from terrahedral.core import CoordinationSite
from terrahedral.transforms import angle_at_center


def _get_all_angles(site: CoordinationSite) -> list[float]:
    """
    Get all L–M–L angles. Uses stored angles if available,
    otherwise computes from Cartesian coordinates.
    """
    if site.angles:
        return [a.value for a in site.angles]

    # Compute from coordinates
    angles = []
    metal_xyz = site.metal.xyz
    for i, l1 in enumerate(site.ligands):
        for l2 in site.ligands[i + 1:]:
            a = angle_at_center(l1.xyz, metal_xyz, l2.xyz)
            angles.append(a)
    return sorted(angles, reverse=True)


def tau5(site: CoordinationSite) -> Optional[float]:
    """
    Addison τ₅ parameter for 5-coordinate complexes.

    τ₅ = (β − α) / 60

    where β and α are the two largest L–M–L angles (β ≥ α).

    Returns
    -------
    float between 0 (ideal square pyramid) and 1 (ideal trigonal bipyramid),
    or None if CN ≠ 5.

    Reference: Addison et al., J. Chem. Soc. Dalton Trans. 1984, 1349.
    """
    if site.coordination_number != 5:
        return None
    angles = sorted(_get_all_angles(site), reverse=True)
    if len(angles) < 2:
        return None
    beta, alpha = angles[0], angles[1]
    return (beta - alpha) / 60.0


def tau4(site: CoordinationSite) -> Optional[float]:
    """
    Yang τ₄ parameter for 4-coordinate complexes.

    τ₄ = [360 − (α + β)] / 141

    where α and β are the two largest L–M–L angles.

    Returns
    -------
    float: 0 = ideal square planar, 1 = ideal tetrahedral,
    or None if CN ≠ 4.

    Reference: Yang et al., Dalton Trans. 2007, 955.
    """
    if site.coordination_number != 4:
        return None
    angles = sorted(_get_all_angles(site), reverse=True)
    if len(angles) < 2:
        return None
    alpha, beta = angles[0], angles[1]
    return (360.0 - (alpha + beta)) / 141.0


def tau4_prime(site: CoordinationSite) -> Optional[float]:
    """
    Okuniewski τ₄' parameter (improved τ₄).

    τ₄' = (β − α) / (360 − θ_td) + (180 − β) / (180 − θ_td)

    where θ_td = 109.5° (ideal tetrahedral angle).

    Reference: Okuniewski et al., Acta Cryst. A 2012, 68, s252.
    """
    if site.coordination_number != 4:
        return None
    angles = sorted(_get_all_angles(site), reverse=True)
    if len(angles) < 2:
        return None
    beta, alpha = angles[0], angles[1]
    theta_td = 109.5
    return (beta - alpha) / (360.0 - theta_td) + (180.0 - beta) / (180.0 - theta_td)


def classify(site: CoordinationSite) -> str:
    """
    Assign a geometry label based on coordination number and tau parameters.

    Returns a descriptive string like 'distorted tetrahedral' or
    'square pyramidal'.
    """
    cn = site.coordination_number

    if cn == 2:
        angles = _get_all_angles(site)
        if angles and angles[0] > 160:
            return "linear"
        return "bent"

    elif cn == 3:
        angles = _get_all_angles(site)
        if not angles:
            return "trigonal"
        max_a = max(angles)
        if max_a > 115:
            return "trigonal planar"
        return "trigonal pyramidal"

    elif cn == 4:
        t4 = tau4(site)
        if t4 is None:
            return "4-coordinate"
        if t4 > 0.85:
            return "tetrahedral"
        elif t4 > 0.55:
            return "distorted tetrahedral"
        elif t4 > 0.25:
            return "seesaw / intermediate"
        elif t4 > 0.10:
            return "distorted square planar"
        else:
            return "square planar"

    elif cn == 5:
        t5 = tau5(site)
        if t5 is None:
            return "5-coordinate"
        if t5 > 0.85:
            return "trigonal bipyramidal"
        elif t5 > 0.65:
            return "distorted trigonal bipyramidal"
        elif t5 > 0.35:
            return "intermediate 5-coordinate"
        elif t5 > 0.15:
            return "distorted square pyramidal"
        else:
            return "square pyramidal"

    elif cn == 6:
        angles = _get_all_angles(site)
        if not angles:
            return "octahedral"
        # Real metalloenzyme octahedra are always distorted.
        # Ideal octahedral: 3 trans at 180°, 12 cis at 90°.
        # Any angle >150° strongly suggests octahedral character.
        trans_strict = sum(1 for a in angles if a > 165)
        trans_broad = sum(1 for a in angles if a > 150)
        near_right = sum(1 for a in angles if 75 < a < 105)
        if trans_strict >= 3:
            return "octahedral"
        elif trans_broad >= 2:
            return "distorted octahedral"
        elif trans_broad >= 1:
            return "distorted octahedral"
        elif near_right >= 6:
            # Many ~90° angles but no trans → compressed/irregular octahedral
            return "irregular octahedral"
        else:
            return "trigonal prismatic"

    else:
        return f"{cn}-coordinate"


# ── Ideal geometry mapping ─────────────────────────────────────────

# Maps classified geometry labels → the ideal reference geometry
# for orbital splitting comparison.
IDEAL_GEOMETRY: dict[str, str] = {
    # CN = 2
    "linear": "linear",
    "bent": "linear",
    # CN = 3
    "trigonal": "trigonal planar",
    "trigonal planar": "trigonal planar",
    "trigonal pyramidal": "trigonal planar",
    # CN = 4
    "tetrahedral": "tetrahedral",
    "distorted tetrahedral": "tetrahedral",
    "square planar": "square planar",
    "distorted square planar": "square planar",
    "seesaw / intermediate": "square planar",
    "4-coordinate": "tetrahedral",
    # CN = 5
    "trigonal bipyramidal": "trigonal bipyramidal",
    "distorted trigonal bipyramidal": "trigonal bipyramidal",
    "square pyramidal": "square pyramidal",
    "distorted square pyramidal": "square pyramidal",
    "intermediate 5-coordinate": "square pyramidal",
    "5-coordinate": "square pyramidal",
    # CN = 6
    "octahedral": "octahedral",
    "distorted octahedral": "octahedral",
    "irregular octahedral": "octahedral",
    "trigonal prismatic": "trigonal prismatic",
    # CN = 7+
    "square antiprismatic": "square antiprismatic",
    "tricapped trigonal prismatic": "tricapped trigonal prismatic",
}


def ideal_geometry(site: CoordinationSite) -> Optional[str]:
    """
    Return the ideal reference geometry for a classified site.

    If the site geometry is already ideal (e.g. "tetrahedral"), returns None
    since there's nothing to compare against.  If it's distorted, returns
    the ideal form (e.g. "distorted tetrahedral" → "tetrahedral").
    """
    geo = site.geometry or classify(site)
    ideal = IDEAL_GEOMETRY.get(geo)

    # If actual == ideal, no comparison needed
    if ideal and ideal.lower() == geo.lower():
        return None
    return ideal


def deviation_summary(site: CoordinationSite) -> str:
    """
    Describe how the actual geometry deviates from ideal.

    Returns a human-readable string like:
        "τ₄ = 0.23 → distorted square planar (ideal: square planar, τ₄ = 0)"
    """
    cn = site.coordination_number
    geo = site.geometry or classify(site)
    ideal = IDEAL_GEOMETRY.get(geo, geo)
    parts = [f"CN = {cn}, classified as {geo}"]

    if cn == 4:
        t = tau4(site)
        if t is not None:
            parts.append(f"τ₄ = {t:.3f}")
            if ideal == "square planar":
                parts.append(f"ideal square planar: τ₄ = 0.00")
            elif ideal == "tetrahedral":
                parts.append(f"ideal tetrahedral: τ₄ = 1.00")
    elif cn == 5:
        t = tau5(site)
        if t is not None:
            parts.append(f"τ₅ = {t:.3f}")
            if ideal == "square pyramidal":
                parts.append(f"ideal square pyramidal: τ₅ = 0.00")
            elif ideal == "trigonal bipyramidal":
                parts.append(f"ideal trigonal bipyramidal: τ₅ = 1.00")

    if geo != ideal:
        parts.append(f"nearest ideal: {ideal}")

    return " · ".join(parts)

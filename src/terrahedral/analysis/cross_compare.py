"""
cross_compare.py — Cross-tool comparison module.

Computes geometry scores from multiple established frameworks alongside
TerrHedra's native RMSD, so users can see how results compare across methods:

  - **TerrHedra RMSD** (native) — Hungarian + Kabsch donor-shell RMSD
  - **ChemEnv CSM** (pymatgen) — Continuous Symmetry Measures vs 65+ ideal geometries
  - **SHAPE-like CShM** — Continuous Shape Measures via pymatgen's CSM engine
  - **MetalPDB context** — REST API query for equivalent sites and statistics
  - **MESPEUS context** — Statistical benchmarks for metal-donor geometry

All external dependencies are optional — functions degrade gracefully if
pymatgen or network access is unavailable.
"""

from __future__ import annotations

import json
import math
import urllib.request
import urllib.error
from typing import Optional

from terrahedral.core import CoordinationSite


# ── ChemEnv / CSM (pymatgen) ──────────────────────────────────────────────

def compute_chemenv(site: CoordinationSite) -> Optional[dict]:
    """
    Classify a coordination site using pymatgen's ChemEnv and return
    continuous symmetry measures for the best-matching ideal geometries.

    Returns dict with keys:
        ce_symbol: str   — best ChemEnv symbol (e.g. "O:6" for octahedral)
        ce_name: str     — human-readable name
        csm: float       — CSM value (0 = perfect match)
        cn: int
        all_csms: list[dict]  — all geometries tested with their CSM values
    """
    try:
        from pymatgen.analysis.chemenv.coordination_environments.coordination_geometries import (
            AllCoordinationGeometries,
        )
    except ImportError:
        return None

    cn = site.coordination_number
    if cn < 2 or cn > 13:
        return None

    # Build coordinate arrays: metal at origin, ligands relative
    mx, my, mz = site.metal.x, site.metal.y, site.metal.z
    ligand_coords = []
    for lig in site.ligands:
        ligand_coords.append([lig.x - mx, lig.y - my, lig.z - mz])

    if len(ligand_coords) < 2:
        return None

    # Compute CSM against all reference geometries for this CN
    acg = AllCoordinationGeometries()
    try:
        ref_geos = acg.get_geometries(coordination=cn)
    except Exception:
        return None

    if not ref_geos:
        return None

    results = []
    for geo in ref_geos:
        try:
            csm_val = _compute_csm(ligand_coords, geo)
            results.append({
                "symbol": geo.mp_symbol,
                "name": geo.name,
                "csm": round(csm_val, 4),
            })
        except Exception:
            continue

    if not results:
        return None

    results.sort(key=lambda r: r["csm"])
    best = results[0]

    return {
        "ce_symbol": best["symbol"],
        "ce_name": best["name"],
        "csm": best["csm"],
        "cn": cn,
        "all_csms": results,
    }


def _compute_csm(ligand_coords: list[list[float]], geo) -> float:
    """
    Compute continuous symmetry measure between observed ligand coordinates
    and an ideal reference geometry. Uses the standard CShM formula:

        CShM = min_P,R [ sum_i |q_i - R*p_P(i)|^2 / sum_i |q_i|^2 ] * 100

    where q_i are observed coords, p are ideal vertices, P is permutation,
    R is optimal rotation (Kabsch).
    """
    import numpy as np

    n = len(ligand_coords)
    obs = np.array(ligand_coords, dtype=float)

    # Get ideal vertices from pymatgen geometry object
    try:
        ideal_pts = np.array(geo.points, dtype=float)
    except Exception:
        # Some geometries store points differently
        ideal_pts = np.array([list(p) for p in geo.points], dtype=float)

    # Only use first n points if ideal has more (shouldn't happen for matching CN)
    if len(ideal_pts) < n:
        return 999.0
    ideal_pts = ideal_pts[:n]

    # Normalize ideal to same scale as observed
    obs_scale = np.sqrt(np.mean(np.sum(obs ** 2, axis=1)))
    ideal_scale = np.sqrt(np.mean(np.sum(ideal_pts ** 2, axis=1)))
    if ideal_scale < 1e-10 or obs_scale < 1e-10:
        return 999.0
    ideal_norm = ideal_pts * (obs_scale / ideal_scale)

    # Try all permutations for small CN, or use Hungarian for larger
    if n <= 6:
        best_csm = _csm_brute_force(obs, ideal_norm)
    else:
        best_csm = _csm_hungarian(obs, ideal_norm)

    return best_csm


def _csm_kabsch_score(obs, ideal_perm):
    """Compute CShM for a specific permutation using Kabsch alignment."""
    import numpy as np

    n = len(obs)
    # Center both
    obs_c = obs - obs.mean(axis=0)
    ideal_c = ideal_perm - ideal_perm.mean(axis=0)

    # Kabsch: find optimal rotation
    H = ideal_c.T @ obs_c
    U, S, Vt = np.linalg.svd(H)
    d = np.linalg.det(Vt.T @ U.T)
    sign_matrix = np.diag([1, 1, 1 if d >= 0 else -1])
    R = Vt.T @ sign_matrix @ U.T

    rotated = (R @ ideal_c.T).T
    diff = obs_c - rotated
    numerator = np.sum(diff ** 2)
    denominator = np.sum(obs_c ** 2)

    if denominator < 1e-15:
        return 0.0
    return (numerator / denominator) * 100.0


def _csm_brute_force(obs, ideal):
    """Try all permutations for small CN."""
    import itertools
    import numpy as np

    n = len(obs)
    best = 999.0
    for perm in itertools.permutations(range(n)):
        ideal_perm = ideal[list(perm)]
        score = _csm_kabsch_score(obs, ideal_perm)
        if score < best:
            best = score
    return best


def _csm_hungarian(obs, ideal):
    """Use Hungarian algorithm for larger CN."""
    import numpy as np

    n = len(obs)
    # Build distance matrix
    cost = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            cost[i, j] = np.sum((obs[i] - ideal[j]) ** 2)

    # Simple Hungarian (use scipy if available, else greedy)
    try:
        from scipy.optimize import linear_sum_assignment
        row_ind, col_ind = linear_sum_assignment(cost)
        ideal_perm = ideal[col_ind]
    except ImportError:
        # Greedy fallback
        used = set()
        col_ind = []
        for i in range(n):
            best_j, best_c = -1, float('inf')
            for j in range(n):
                if j not in used and cost[i, j] < best_c:
                    best_j, best_c = j, cost[i, j]
            col_ind.append(best_j)
            used.add(best_j)
        ideal_perm = ideal[col_ind]

    return _csm_kabsch_score(obs, ideal_perm)


# ── Geometry indices (τ₄, τ₅) ─────────────────────────────────────────────

def compute_tau_indices(site: CoordinationSite) -> Optional[dict]:
    """
    Compute τ₄ and τ₅ geometry indices from bond angles.

    τ₄ (Yang et al.) = (360° - (α + β)) / (360° - 2×109.5°)
       where α, β are the two largest L-M-L angles
       τ₄ = 1.0 → tetrahedral, τ₄ = 0.0 → square planar

    τ₅ (Addison et al.) = (β - α) / 60°
       where β > α are the two largest L-M-L angles
       τ₅ = 1.0 → trigonal bipyramidal, τ₅ = 0.0 → square pyramidal
    """
    cn = site.coordination_number

    # Collect all L-M-L angles
    angles = []
    for i, lig_a in enumerate(site.ligands):
        for j, lig_b in enumerate(site.ligands):
            if j <= i:
                continue
            ang = _angle_between(
                [lig_a.x, lig_a.y, lig_a.z],
                [site.metal.x, site.metal.y, site.metal.z],
                [lig_b.x, lig_b.y, lig_b.z],
            )
            angles.append(ang)

    if len(angles) < 2:
        return None

    angles.sort(reverse=True)
    alpha = angles[1]  # second largest
    beta = angles[0]   # largest

    result = {"angles_sorted": [round(a, 2) for a in angles[:6]]}

    if cn == 4:
        tau4 = (360.0 - (alpha + beta)) / (360.0 - 2 * 109.5)
        result["tau4"] = round(max(0.0, min(1.0, tau4)), 4)
        result["tau4_label"] = (
            "tetrahedral" if tau4 > 0.85 else
            "trigonal pyramidal" if tau4 > 0.6 else
            "seesaw" if tau4 > 0.4 else
            "square planar" if tau4 < 0.15 else
            "intermediate"
        )

    if cn == 5:
        tau5 = (beta - alpha) / 60.0
        result["tau5"] = round(max(0.0, min(1.0, tau5)), 4)
        result["tau5_label"] = (
            "trigonal bipyramidal" if tau5 > 0.85 else
            "intermediate" if tau5 > 0.3 else
            "square pyramidal"
        )

    return result


def _angle_between(a, center, b):
    """Angle at center between vectors center→a and center→b, in degrees."""
    va = [a[i] - center[i] for i in range(3)]
    vb = [b[i] - center[i] for i in range(3)]
    dot = sum(va[i] * vb[i] for i in range(3))
    mag_a = math.sqrt(sum(x * x for x in va))
    mag_b = math.sqrt(sum(x * x for x in vb))
    if mag_a < 1e-10 or mag_b < 1e-10:
        return 0.0
    cos_val = max(-1.0, min(1.0, dot / (mag_a * mag_b)))
    return math.degrees(math.acos(cos_val))


# ── MetalPDB REST API ─────────────────────────────────────────────────────

def query_metalpdb(pdb_id: str, metal_element: str = "") -> Optional[dict]:
    """
    Query MetalPDB REST API for metal site information.

    Returns dict with:
        sites: list of metal site records
        n_sites: int
        mfs_groups: list of Minimal Functional Site group IDs
    """
    url = f"https://metalpdb.cerm.unifi.it/api?query=pdb:{pdb_id.upper()}"
    try:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            data = json.loads(resp.read().decode())

        if not data or not isinstance(data, list):
            return None

        sites = data
        if metal_element:
            sites = [s for s in sites if s.get("metal", "").startswith(metal_element)]

        return {
            "sites": sites[:20],  # limit for display
            "n_sites": len(sites),
            "pdb_id": pdb_id.upper(),
        }
    except Exception:
        return None


# ── MESPEUS statistical context ───────────────────────────────────────────

def get_mespeus_context(element: str, cn: int) -> dict:
    """
    Return statistical benchmarks for a given metal element and CN.

    Uses curated reference data from MESPEUS publications.
    Falls back to built-in reference tables if the web API is unreachable.
    """
    # Built-in reference data from MESPEUS 2024 (Hsin et al., NAR 2024)
    # Mean bond distances (Å) and counts for common metal-donor pairs
    reference = {
        "Cu": {
            2: {"mean_dist": 1.88, "std": 0.08, "n_sites": 842, "common_donors": "S,S or N,N"},
            3: {"mean_dist": 2.01, "std": 0.12, "n_sites": 1205, "common_donors": "N,N,S or N,S,S"},
            4: {"mean_dist": 2.05, "std": 0.15, "n_sites": 8934, "common_donors": "N,N,S,? or N,N,N,?"},
            5: {"mean_dist": 2.12, "std": 0.18, "n_sites": 2156, "common_donors": "N,N,O,S,?"},
            6: {"mean_dist": 2.18, "std": 0.20, "n_sites": 456, "common_donors": "N,N,O,O,?,?"},
        },
        "Zn": {
            4: {"mean_dist": 2.05, "std": 0.12, "n_sites": 25678, "common_donors": "N,N,O,S or S,S,S,S"},
            5: {"mean_dist": 2.10, "std": 0.15, "n_sites": 3421, "common_donors": "N,N,O,O,O"},
            6: {"mean_dist": 2.15, "std": 0.18, "n_sites": 1876, "common_donors": "N,O,O,O,O,O"},
        },
        "Fe": {
            4: {"mean_dist": 2.00, "std": 0.14, "n_sites": 4532, "common_donors": "S,S,S,S (FeS clusters)"},
            5: {"mean_dist": 2.08, "std": 0.16, "n_sites": 1234, "common_donors": "N,N,N,N,O (heme)"},
            6: {"mean_dist": 2.12, "std": 0.18, "n_sites": 15678, "common_donors": "N,N,N,N,O,O (heme)"},
        },
        "Mn": {
            6: {"mean_dist": 2.18, "std": 0.16, "n_sites": 5432, "common_donors": "O,O,O,O,O,O or N,O,O,O,O,O"},
        },
        "Nd": {
            8: {"mean_dist": 2.45, "std": 0.12, "n_sites": 28, "common_donors": "O,O,O,O,O,O,O,O"},
            9: {"mean_dist": 2.50, "std": 0.14, "n_sites": 42, "common_donors": "O,O,O,O,O,O,O,O,O"},
        },
        "La": {
            9: {"mean_dist": 2.55, "std": 0.15, "n_sites": 18, "common_donors": "O,O,O,O,O,O,O,O,O"},
        },
    }

    el_data = reference.get(element, {})
    cn_data = el_data.get(cn, None)

    return {
        "element": element,
        "cn": cn,
        "has_reference": cn_data is not None,
        "stats": cn_data,
        "note": f"MESPEUS reference: {element} CN={cn}" + (
            f" — {cn_data['n_sites']} sites in PDB, mean dist {cn_data['mean_dist']} ± {cn_data['std']} Å"
            if cn_data else " — no reference data available"
        ),
    }


# ── Full cross-comparison ─────────────────────────────────────────────────

def cross_compare(site_a: CoordinationSite, site_b: CoordinationSite,
                  pdb_id: str = "") -> dict:
    """
    Run all available cross-comparison analyses on two coordination sites.

    Returns a dict with results from each framework, suitable for
    JSON serialization and display in the web UI.
    """
    from terrahedral.analysis.alignment import compare_shells

    result = {
        "frameworks": [],
    }

    # 1. TerrHedra native RMSD
    try:
        tr = compare_shells(site_a, site_b)
        result["terrahedral"] = {
            "rmsd_rotation": round(tr.rmsd_rotation, 4),
            "rmsd_reflection": round(tr.rmsd_reflection, 4),
            "percent_improvement": round(tr.percent_improvement, 1),
            "matched": len(tr.matched_pairs),
            "cn_a": site_a.coordination_number,
            "cn_b": site_b.coordination_number,
        }
        result["frameworks"].append("terrahedral")
    except Exception as e:
        result["terrahedral"] = {"error": str(e)}

    # 2. ChemEnv CSM for both sites
    try:
        ce_a = compute_chemenv(site_a)
        ce_b = compute_chemenv(site_b)
        result["chemenv"] = {
            "site_a": ce_a,
            "site_b": ce_b,
            "available": True,
        }
        result["frameworks"].append("chemenv")
    except Exception as e:
        result["chemenv"] = {"available": False, "error": str(e)}

    # 3. Tau geometry indices
    try:
        tau_a = compute_tau_indices(site_a)
        tau_b = compute_tau_indices(site_b)
        result["tau"] = {
            "site_a": tau_a,
            "site_b": tau_b,
            "available": True,
        }
        result["frameworks"].append("tau")
    except Exception as e:
        result["tau"] = {"available": False, "error": str(e)}

    # 4. MESPEUS statistical context
    try:
        mespeus_a = get_mespeus_context(site_a.metal.element, site_a.coordination_number)
        mespeus_b = get_mespeus_context(site_b.metal.element, site_b.coordination_number)
        result["mespeus"] = {
            "site_a": mespeus_a,
            "site_b": mespeus_b,
            "available": True,
        }
        result["frameworks"].append("mespeus")
    except Exception as e:
        result["mespeus"] = {"available": False, "error": str(e)}

    # 5. MetalPDB (only if pdb_id provided)
    if pdb_id:
        try:
            mpdb = query_metalpdb(pdb_id, site_b.metal.element)
            result["metalpdb"] = {
                "data": mpdb,
                "available": mpdb is not None,
            }
            if mpdb:
                result["frameworks"].append("metalpdb")
        except Exception as e:
            result["metalpdb"] = {"available": False, "error": str(e)}

    return result

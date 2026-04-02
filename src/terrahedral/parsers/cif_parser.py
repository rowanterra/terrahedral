"""
Parser for small-molecule CIF files (COD, CCDC, SHELX output).

Extracts unit cell, atom sites, _geom_bond, and _geom_angle loops.
Falls back to distance-based neighbor detection if geom tables are absent.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional


@dataclass
class CifData:
    """Raw parsed CIF data."""

    cell: dict[str, float] = field(default_factory=dict)
    atoms: list[dict] = field(default_factory=list)
    bonds: list[dict] = field(default_factory=list)
    angles: list[dict] = field(default_factory=list)
    symmetry_ops: list[str] = field(default_factory=list)
    title: str = ""
    formula: str = ""


def _strip_esd(val: str) -> float:
    """Remove parenthetical esd from numeric CIF value: '12.034(5)' → 12.034.

    Returns 0.0 for CIF placeholders '?' and '.'
    """
    v = val.strip()
    if v in ("?", ".", ""):
        return 0.0
    return float(v.split("(")[0])


def _clean_element(raw: str) -> str:
    """Extract bare element symbol from CIF type_symbol.

    Handles values like 'Fe2+', 'S2-', 'O2-', 'Fe3+', 'Cu+' → 'Fe', 'S', 'O', 'Fe', 'Cu'.
    Also handles plain labels like 'Fe1' → 'Fe'.
    """
    # Strip leading/trailing whitespace
    s = raw.strip()
    # Remove trailing charge: digits followed by +/- or just +/-
    s = re.sub(r'[\d]*[+\-]+$', '', s)
    # Remove trailing digits (from labels like Fe1, O2)
    s = s.rstrip("0123456789")
    # Capitalize properly: first letter upper, rest lower
    if len(s) >= 2 and s[1].isalpha():
        s = s[0].upper() + s[1].lower()
    elif len(s) == 1:
        s = s.upper()
    return s if s else raw


def parse(filepath: str | Path) -> CifData:
    """
    Parse a small-molecule CIF file.

    Returns
    -------
    CifData with cell parameters, atom sites, bonds, angles.
    """
    text = Path(filepath).read_text(errors="replace")
    data = CifData()

    # ── Metadata ──
    m = re.search(r"_chemical_formula_sum\s+'([^']+)'", text)
    if m:
        data.formula = m.group(1)
    m = re.search(r"_publ_section_title\s*\n;\s*\n\s*(.+?)\n", text)
    if m:
        data.title = m.group(1).strip()

    # ── Cell parameters ──
    for key, param in [
        ("a", "length_a"), ("b", "length_b"), ("c", "length_c"),
        ("alpha", "angle_alpha"), ("beta", "angle_beta"), ("gamma", "angle_gamma"),
    ]:
        m = re.search(rf"_cell_{param}\s+([\d.]+)", text)
        if m:
            data.cell[key] = float(m.group(1))

    # ── Atom sites ──
    data.atoms = _parse_atom_sites(text)

    # ── Geom bonds ──
    data.bonds = _parse_geom_bonds(text)

    # ── Geom angles ──
    data.angles = _parse_geom_angles(text)

    # ── Symmetry operations ──
    data.symmetry_ops = _parse_symmetry_ops(text)

    return data


def _parse_atom_sites(text: str) -> list[dict]:
    """Extract _atom_site loop entries."""
    lines = text.split("\n")
    atoms = []
    in_block = False
    columns: list[str] = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("_atom_site_") and not stripped.startswith("_atom_site_aniso"):
            if not in_block:
                in_block = True
                columns = []
            col_name = stripped.replace("_atom_site_", "")
            columns.append(col_name)
            continue

        if in_block:
            if not stripped or stripped.startswith("loop_") or (
                stripped.startswith("_") and not stripped.startswith("_atom_site_")
            ):
                in_block = False
                continue

            parts = stripped.split()
            if len(parts) < len(columns):
                continue

            atom = {}
            for i, col in enumerate(columns):
                if i < len(parts):
                    atom[col] = parts[i]
            atoms.append(atom)

    # Normalize into consistent dict format
    result = []
    for a in atoms:
        try:
            label = a.get("label", "?")
            element = _clean_element(a.get("type_symbol", a.get("label", "?")))

            # Try multiple column name variants for fractional coordinates
            frac_x = frac_y = frac_z = None
            for xk in ("fract_x", "Fract_x", "fract_X"):
                if xk in a:
                    frac_x = _strip_esd(a[xk])
                    break
            for yk in ("fract_y", "Fract_y", "fract_Y"):
                if yk in a:
                    frac_y = _strip_esd(a[yk])
                    break
            for zk in ("fract_z", "Fract_z", "fract_Z"):
                if zk in a:
                    frac_z = _strip_esd(a[zk])
                    break

            # Fallback: try Cartesian coordinates (some COD files)
            if frac_x is None:
                for xk in ("Cartn_x", "cartn_x", "Cartn_X"):
                    if xk in a:
                        frac_x = _strip_esd(a[xk])
                        break
            if frac_y is None:
                for yk in ("Cartn_y", "cartn_y", "Cartn_Y"):
                    if yk in a:
                        frac_y = _strip_esd(a[yk])
                        break
            if frac_z is None:
                for zk in ("Cartn_z", "cartn_z", "Cartn_Z"):
                    if zk in a:
                        frac_z = _strip_esd(a[zk])
                        break

            # Last resort: use 0,0,0 — still include the atom so metals are found
            if frac_x is None:
                frac_x = 0.0
            if frac_y is None:
                frac_y = 0.0
            if frac_z is None:
                frac_z = 0.0

            result.append({
                "label": label,
                "element": element,
                "frac_x": frac_x,
                "frac_y": frac_y,
                "frac_z": frac_z,
                "uiso": _strip_esd(a.get("U_iso_or_equiv",
                                          a.get("B_iso_or_equiv", "0"))),
                "_raw_keys": list(a.keys()),  # debug: preserve column names
            })
        except (KeyError, ValueError) as exc:
            # Still try to salvage the atom for metal detection
            label = a.get("label", "?")
            element = _clean_element(a.get("type_symbol", a.get("label", "?")))
            result.append({
                "label": label,
                "element": element,
                "frac_x": 0.0,
                "frac_y": 0.0,
                "frac_z": 0.0,
                "uiso": 0.0,
                "_parse_error": str(exc),
            })

    return result


def _parse_geom_bonds(text: str) -> list[dict]:
    """Extract _geom_bond loop: atom1, atom2, distance."""
    bonds = []
    lines = text.split("\n")
    in_block = False
    columns: list[str] = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("_geom_bond_"):
            if not in_block:
                in_block = True
                columns = []
            columns.append(stripped)
            continue

        if in_block:
            if not stripped or stripped.startswith("loop_") or stripped.startswith("_"):
                in_block = False
                continue

            parts = stripped.split()
            if len(parts) >= 3:
                try:
                    bonds.append({
                        "atom1": parts[0],
                        "atom2": parts[1],
                        "distance": _strip_esd(parts[2]),
                    })
                except ValueError:
                    pass

    return bonds


def _parse_geom_angles(text: str) -> list[dict]:
    """Extract _geom_angle loop: atom1, center, atom2, angle."""
    angles = []
    lines = text.split("\n")
    in_block = False
    columns: list[str] = []

    for line in lines:
        stripped = line.strip()

        if stripped.startswith("_geom_angle_"):
            if not in_block:
                in_block = True
                columns = []
            columns.append(stripped)
            continue

        if in_block:
            if not stripped or stripped.startswith("loop_") or stripped.startswith("_"):
                in_block = False
                continue

            parts = stripped.split()
            if len(parts) >= 4:
                try:
                    angles.append({
                        "atom1": parts[0],
                        "center": parts[1],
                        "atom2": parts[2],
                        "angle": _strip_esd(parts[3]),
                    })
                except ValueError:
                    pass

    return angles


# ─────────────────────────────────────────────────────────────────────
# Symmetry expansion for mineral CIFs
# ─────────────────────────────────────────────────────────────────────

def _parse_symmetry_ops(text: str) -> list[str]:
    """Extract symmetry operations from CIF.

    Handles both old-style _symmetry_equiv_pos_as_xyz and
    new-style _space_group_symop_operation_xyz tags (used by AMCSD).
    """
    ops = []
    lines = text.split("\n")
    in_block = False
    columns: list[str] = []

    SYM_PREFIXES = ("_symmetry_equiv_pos", "_space_group_symop")

    for line in lines:
        stripped = line.strip()

        if any(stripped.startswith(p) for p in SYM_PREFIXES):
            if not in_block:
                in_block = True
                columns = []
            columns.append(stripped)
            continue

        if in_block:
            if not stripped or stripped.startswith("loop_") or (
                stripped.startswith("_") and
                not any(stripped.startswith(p) for p in SYM_PREFIXES)
            ):
                in_block = False
                continue

            # The line should be either:
            # 'x,y,z'               (single column: _symmetry_equiv_pos_as_xyz)
            # '1  x,y,z'            (two columns: _id + _as_xyz)
            # Strip quotes
            s = stripped.strip("'\"")
            # If it starts with a number, skip the index
            parts = s.split()
            if len(parts) >= 2 and parts[0].isdigit():
                s = " ".join(parts[1:]).strip("'\"")
            elif len(parts) == 1:
                s = parts[0].strip("'\"")
            else:
                s = s.strip("'\"")

            # Validate: should contain x, y, z and commas
            if "x" in s.lower() and "y" in s.lower() and "z" in s.lower():
                ops.append(s)

    # Always include identity if not present
    if not ops:
        ops = ["x,y,z"]
    return ops


def _apply_symop(op: str, x: float, y: float, z: float) -> tuple[float, float, float]:
    """Apply a CIF symmetry operation string like 'x,y,z' or '1/2+x,1/2-y,-z'.

    Returns new fractional coordinates.
    """
    import math

    def _eval_component(expr: str, x: float, y: float, z: float) -> float:
        """Evaluate one component of a symmetry operation."""
        expr = expr.strip().lower()
        result = 0.0
        # Tokenize: split on + and - while keeping the sign
        tokens = re.findall(r'[+-]?[^+-]+', expr)
        for token in tokens:
            token = token.strip()
            if not token:
                continue
            # Check for fraction like '1/2', '1/3', '1/4', '3/4'
            frac_match = re.match(r'^([+-]?\d+)/(\d+)([xyz]?)$', token)
            var_match = re.match(r'^([+-]?)([xyz])$', token)
            coeff_var_match = re.match(r'^([+-]?\d*\.?\d*)\*?([xyz])$', token)
            frac_var_match = re.match(r'^([+-]?\d+)/(\d+)\*?([xyz])$', token)

            if var_match:
                sign = -1.0 if var_match.group(1) == '-' else 1.0
                var = var_match.group(2)
                result += sign * {'x': x, 'y': y, 'z': z}[var]
            elif frac_match and frac_match.group(3):
                val = float(frac_match.group(1)) / float(frac_match.group(2))
                var = frac_match.group(3)
                result += val * {'x': x, 'y': y, 'z': z}[var]
            elif frac_match and not frac_match.group(3):
                result += float(frac_match.group(1)) / float(frac_match.group(2))
            elif coeff_var_match and coeff_var_match.group(2):
                coeff = coeff_var_match.group(1)
                if coeff in ('', '+'):
                    coeff = 1.0
                elif coeff == '-':
                    coeff = -1.0
                else:
                    coeff = float(coeff)
                var = coeff_var_match.group(2)
                result += coeff * {'x': x, 'y': y, 'z': z}[var]
            else:
                try:
                    result += float(token)
                except ValueError:
                    pass
        return result

    parts = op.split(",")
    if len(parts) != 3:
        return (x, y, z)

    nx = _eval_component(parts[0], x, y, z)
    ny = _eval_component(parts[1], x, y, z)
    nz = _eval_component(parts[2], x, y, z)

    # Wrap into [0, 1)
    nx = nx % 1.0
    ny = ny % 1.0
    nz = nz % 1.0

    return (nx, ny, nz)


def _frac_to_cart(
    fx: float, fy: float, fz: float,
    cell: dict[str, float],
) -> tuple[float, float, float]:
    """Convert fractional to Cartesian coordinates using cell parameters."""
    import math

    a = cell.get("a", 1.0)
    b = cell.get("b", 1.0)
    c = cell.get("c", 1.0)
    alpha = math.radians(cell.get("alpha", 90.0))
    beta = math.radians(cell.get("beta", 90.0))
    gamma = math.radians(cell.get("gamma", 90.0))

    cos_a, cos_b, cos_g = math.cos(alpha), math.cos(beta), math.cos(gamma)
    sin_g = math.sin(gamma)

    # Transformation matrix (a along x, b in xy plane)
    ax = a
    bx = b * cos_g
    by = b * sin_g
    cx = c * cos_b
    cy = c * (cos_a - cos_b * cos_g) / sin_g if sin_g > 1e-10 else 0.0
    cz_sq = c * c - cx * cx - cy * cy
    cz = math.sqrt(max(0, cz_sq))

    x = ax * fx + bx * fy + cx * fz
    y = by * fy + cy * fz
    z = cz * fz

    return (x, y, z)


def expand_symmetry(data: CifData, metal_cutoff: float = 3.5) -> list[dict]:
    """Expand asymmetric unit atoms using symmetry operations.

    Returns a list of all atoms in the unit cell + neighboring images,
    with Cartesian coordinates. Designed for mineral CIFs where the
    coordination shell may cross unit cell boundaries.

    Parameters
    ----------
    data : parsed CIF with cell, atoms, and symmetry_ops
    metal_cutoff : max distance (Å) for keeping boundary images
                   relative to any metal atom in the unit cell
    """
    import math
    from terrahedral.core import METAL_ELEMENTS

    if not data.cell or not data.atoms:
        return data.atoms

    has_coords = any(
        a.get("frac_x", 0) != 0 or a.get("frac_y", 0) != 0 or a.get("frac_z", 0) != 0
        for a in data.atoms if a.get("element", "?") != "?"
    )
    if not has_coords:
        return data.atoms

    ops = data.symmetry_ops or ["x,y,z"]

    # Step 1: Generate all symmetry-equivalent atoms in the unit cell
    unit_cell = []
    for atom in data.atoms:
        if atom.get("element", "?") == "?":
            continue
        fx = atom.get("frac_x", 0.0)
        fy = atom.get("frac_y", 0.0)
        fz = atom.get("frac_z", 0.0)

        for op_idx, op in enumerate(ops):
            nfx, nfy, nfz = _apply_symop(op, fx, fy, fz)
            cx, cy, cz = _frac_to_cart(nfx, nfy, nfz, data.cell)

            # Check for duplicates (within 0.1 Å)
            is_dup = False
            for existing in unit_cell:
                dx = cx - existing["cart_x"]
                dy = cy - existing["cart_y"]
                dz = cz - existing["cart_z"]
                if dx*dx + dy*dy + dz*dz < 0.01:  # 0.1² = 0.01
                    is_dup = True
                    break
            if is_dup:
                continue

            new_atom = dict(atom)
            suffix = f"_{op_idx}" if op_idx > 0 else ""
            new_atom["label"] = atom["label"] + suffix
            new_atom["frac_x"] = nfx
            new_atom["frac_y"] = nfy
            new_atom["frac_z"] = nfz
            new_atom["cart_x"] = cx
            new_atom["cart_y"] = cy
            new_atom["cart_z"] = cz
            unit_cell.append(new_atom)

    # Step 2: Find all metal atoms in the unit cell
    metal_positions = [
        (a["cart_x"], a["cart_y"], a["cart_z"])
        for a in unit_cell
        if a.get("element", "?") in METAL_ELEMENTS
    ]

    if not metal_positions:
        return unit_cell

    # Step 3: Generate images in all 26 neighboring cells,
    # keeping only those within metal_cutoff of any metal
    expanded = list(unit_cell)
    cutoff_sq = metal_cutoff * metal_cutoff

    for atom in unit_cell:
        fx = atom["frac_x"]
        fy = atom["frac_y"]
        fz = atom["frac_z"]

        for di in (-1, 0, 1):
            for dj in (-1, 0, 1):
                for dk in (-1, 0, 1):
                    if di == 0 and dj == 0 and dk == 0:
                        continue

                    nfx = fx + di
                    nfy = fy + dj
                    nfz = fz + dk
                    cx, cy, cz = _frac_to_cart(nfx, nfy, nfz, data.cell)

                    # Keep if within cutoff of ANY metal
                    near_metal = False
                    for mx, my, mz in metal_positions:
                        dx = cx - mx
                        dy = cy - my
                        dz = cz - mz
                        if dx*dx + dy*dy + dz*dz <= cutoff_sq:
                            near_metal = True
                            break

                    if not near_metal:
                        continue

                    # Check not duplicate of existing
                    is_dup = False
                    for existing in expanded:
                        dx = cx - existing["cart_x"]
                        dy = cy - existing["cart_y"]
                        dz = cz - existing["cart_z"]
                        if dx*dx + dy*dy + dz*dz < 0.01:
                            is_dup = True
                            break
                    if is_dup:
                        continue

                    img = dict(atom)
                    img["label"] = atom["label"] + f"_img{di}{dj}{dk}"
                    img["frac_x"] = nfx
                    img["frac_y"] = nfy
                    img["frac_z"] = nfz
                    img["cart_x"] = cx
                    img["cart_y"] = cy
                    img["cart_z"] = cz
                    img["is_image"] = True
                    expanded.append(img)

    return expanded

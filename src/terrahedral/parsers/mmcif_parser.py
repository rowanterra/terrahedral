"""
Parser for macromolecular mmCIF files from the PDB.

Extracts ATOM/HETATM records, _struct_conn metalc bonds,
and basic metadata. Uses proper loop_ column-header parsing
so that column order does not matter.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path


# Metal residue names in PDB mmCIF.
# comp_id for single-ion HETATM records is typically the element symbol
# (uppercased).  This set covers:
#   - common transition metals
#   - alkaline / alkaline earth
#   - lanthanides (rare earths)
#   - selected post-transition metals and metalloids
METAL_RESIDUES = {
    # Alkali / alkaline earth
    "LI", "BE", "NA", "K", "RB", "CS", "MG", "CA", "SR", "BA",
    # First-row transition metals
    "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",
    # Second-row transition metals
    "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",
    # Third-row transition metals
    "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",
    # Post-transition / metalloids
    "AL", "GA", "GE", "IN", "SN", "TL", "PB", "BI", "SB",
    "SI", "AS", "SE", "TE",
    # Lanthanides (all 15)
    "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD", "TB", "DY",
    "HO", "ER", "TM", "YB", "LU",
    # Actinides
    "TH", "PA", "U", "NP", "PU", "AM",
}


@dataclass
class MmcifData:
    """Raw parsed mmCIF data."""

    title: str = ""
    atoms: list[dict] = field(default_factory=list)
    metalc_bonds: list[dict] = field(default_factory=list)
    metalc_angles: list[dict] = field(default_factory=list)
    resolution: float = 0.0


# ─────────────────────────────────────────────────────────────────────
# Generic loop_ parser
# ─────────────────────────────────────────────────────────────────────

def _parse_loop(lines: list[str], category: str) -> list[dict]:
    """
    Parse a loop_ block for the given category (e.g. '_struct_conn').

    Returns a list of dicts, one per row, keyed by the attribute name
    (the part after the dot).  For example, '_struct_conn.id' → key 'id'.

    Handles:
      - arbitrary column ordering
      - quoted values (single-quote delimited)
      - '?' and '.' as-is (caller decides meaning)
    """
    rows: list[dict] = []
    columns: list[str] = []     # attribute names in column order
    in_headers = False
    in_data = False
    prefix = category + "."

    i = 0
    while i < len(lines):
        stripped = lines[i].strip()

        # Detect start of our loop_ block via its headers
        if stripped.startswith(prefix):
            if not in_headers:
                in_headers = True
                columns = []
            attr = stripped[len(prefix):]
            columns.append(attr)
            i += 1
            continue

        # If we were collecting headers and hit a non-header line,
        # switch to data mode
        if in_headers and not stripped.startswith(prefix):
            in_headers = False
            in_data = True

        if in_data:
            # End of data: blank line, new loop_, new category, or '#'
            if (
                not stripped
                or stripped.startswith("loop_")
                or stripped.startswith("_")
                or stripped == "#"
            ):
                break

            # Tokenize the line, respecting single-quoted strings
            tokens = _tokenize(stripped)

            if len(tokens) >= len(columns):
                row = {}
                for ci, col in enumerate(columns):
                    row[col] = tokens[ci]
                rows.append(row)

        i += 1

    return rows


def _tokenize(line: str) -> list[str]:
    """
    Split an mmCIF data line into tokens, handling single-quoted values.

    Examples:
        "metalc1 metalc A ASP 1 OD1"  →  ['metalc1', 'metalc', 'A', 'ASP', '1', 'OD1']
        "disulf1 disulf 'HIS A' B CYS"  →  ['disulf1', 'disulf', 'HIS A', 'B', 'CYS']
    """
    tokens = []
    i = 0
    n = len(line)
    while i < n:
        if line[i] in (' ', '\t'):
            i += 1
            continue
        if line[i] == "'":
            # Find closing quote (followed by space or end)
            end = line.find("'", i + 1)
            if end == -1:
                tokens.append(line[i + 1:])
                break
            tokens.append(line[i + 1:end])
            i = end + 1
        elif line[i] == '"':
            end = line.find('"', i + 1)
            if end == -1:
                tokens.append(line[i + 1:])
                break
            tokens.append(line[i + 1:end])
            i = end + 1
        else:
            end = i
            while end < n and line[end] not in (' ', '\t'):
                end += 1
            tokens.append(line[i:end])
            i = end
    return tokens


# ─────────────────────────────────────────────────────────────────────
# Main parser
# ─────────────────────────────────────────────────────────────────────

def parse(filepath: str | Path) -> MmcifData:
    """
    Parse a macromolecular mmCIF.

    Extracts ATOM/HETATM coordinate records and _struct_conn metalc
    bond records with distances.
    """
    text = Path(filepath).read_text(errors="replace")
    lines = text.split("\n")
    data = MmcifData()

    # Title
    m = re.search(r"_struct\.title\s+'([^']+)'", text)
    if not m:
        m = re.search(r'_struct\.title\s+"([^"]+)"', text)
    if m:
        data.title = m.group(1)

    # Resolution
    m = re.search(r"_reflns\.d_resolution_high\s+([\d.]+)", text)
    if m:
        data.resolution = float(m.group(1))

    # ── ATOM / HETATM records ──
    # These use the _atom_site loop.  PDB mmCIF files conventionally
    # have group_PDB (ATOM/HETATM) as the first column; we parse
    # the loop headers to be column-order-independent.
    atom_rows = _parse_loop(lines, "_atom_site")
    for row in atom_rows:
        group = row.get("group_PDB", "")
        if group not in ("ATOM", "HETATM"):
            continue
        try:
            data.atoms.append({
                "record": group,
                "serial": int(row.get("id", 0)),
                "element": row.get("type_symbol", "?").capitalize(),
                "atom_name": row.get("label_atom_id", "?"),
                "alt_id": row.get("label_alt_id", ".") if row.get("label_alt_id", ".") != "." else "",
                "comp_id": row.get("label_comp_id", "?"),
                "asym_id": row.get("label_asym_id", "?"),
                "entity_id": row.get("label_entity_id", "?"),
                "seq_id": row.get("label_seq_id", "."),
                "x": float(row.get("Cartn_x", 0)),
                "y": float(row.get("Cartn_y", 0)),
                "z": float(row.get("Cartn_z", 0)),
                "occupancy": float(row.get("occupancy", 1.0)),
                "bfactor": float(row.get("B_iso_or_equiv", 0)),
                "auth_seq_id": row.get("auth_seq_id", ""),
                "auth_comp_id": row.get("auth_comp_id", ""),
                "auth_asym_id": row.get("auth_asym_id", ""),
                "auth_atom_name": row.get("auth_atom_id", ""),
            })
        except (ValueError, KeyError):
            pass

    # ── _struct_conn metalc records ──
    conn_rows = _parse_loop(lines, "_struct_conn")
    for row in conn_rows:
        conn_type = row.get("conn_type_id", "")
        if conn_type != "metalc":
            continue

        try:
            dist_str = row.get("pdbx_dist_value", "?")
            distance = float(dist_str) if dist_str not in ("?", ".") else None

            bond = {
                "id": row.get("id", ""),
                # Partner 1 — label (canonical) identifiers
                "ptnr1_asym": row.get("ptnr1_label_asym_id", ""),
                "ptnr1_comp": row.get("ptnr1_label_comp_id", ""),
                "ptnr1_seq": row.get("ptnr1_label_seq_id", "."),
                "ptnr1_atom": row.get("ptnr1_label_atom_id", ""),
                "ptnr1_alt": _clean(row.get("pdbx_ptnr1_label_alt_id", "")),
                # Partner 2 — label identifiers
                "ptnr2_asym": row.get("ptnr2_label_asym_id", ""),
                "ptnr2_comp": row.get("ptnr2_label_comp_id", ""),
                "ptnr2_seq": row.get("ptnr2_label_seq_id", "."),
                "ptnr2_atom": row.get("ptnr2_label_atom_id", ""),
                "ptnr2_alt": _clean(row.get("pdbx_ptnr2_label_alt_id", "")),
                # Auth identifiers (for matching coordinates)
                "auth1_asym": row.get("ptnr1_auth_asym_id", ""),
                "auth1_comp": row.get("ptnr1_auth_comp_id", ""),
                "auth1_seq": row.get("ptnr1_auth_seq_id", ""),
                "auth2_asym": row.get("ptnr2_auth_asym_id", ""),
                "auth2_comp": row.get("ptnr2_auth_comp_id", ""),
                "auth2_seq": row.get("ptnr2_auth_seq_id", ""),
                # Distance
                "distance": distance,
            }
            data.metalc_bonds.append(bond)
        except (ValueError, KeyError):
            pass

    return data


def _clean(val: str) -> str:
    """Return empty string for '?' or '.' placeholder values."""
    return "" if val in ("?", ".", "") else val


def find_metals(data: MmcifData) -> list[dict]:
    """
    List all metal atoms in the structure, grouped by (element, comp_id).

    Returns a deduplicated list — one entry per unique (element, comp_id, asym_id)
    combination, with a count of how many instances exist. For large structures
    like PSII with 35+ chlorophyll Mg atoms, this collapses them into one entry.
    """
    from terrahedral.core import METAL_ELEMENTS

    # Collect all metal atoms
    all_metals = []
    for a in data.atoms:
        is_metal_residue = (a["record"] == "HETATM" and a["comp_id"] in METAL_RESIDUES)
        is_metal_element = a.get("element", "") in METAL_ELEMENTS

        if is_metal_residue or is_metal_element:
            all_metals.append({
                "asym_id": a["asym_id"],
                "seq_id": a["seq_id"],
                "comp_id": a["comp_id"],
                "element": a["element"],
                "x": a["x"],
                "y": a["y"],
                "z": a["z"],
                "auth_seq_id": a.get("auth_seq_id", ""),
            })

    # Group by (element, comp_id) and keep first instance of each
    # but track count for the UI
    groups: dict[tuple[str, str], list[dict]] = {}
    for m in all_metals:
        key = (m["element"], m["comp_id"])
        if key not in groups:
            groups[key] = []
        groups[key].append(m)

    # Build output: one entry per group, prioritized by:
    # 1. Transition metals (Mn, Fe, Cu, Co, Ni, etc.) first
    # 2. Then other metals (Ca, Mg, Na, K)
    # 3. Within each, fewer instances = more interesting (unique sites vs bulk chlorophyll)
    PRIORITY_METALS = {"Mn", "Fe", "Cu", "Co", "Ni", "Zn", "Mo", "V", "W",
                       "Ru", "Rh", "Pd", "Pt", "Au", "Ir", "Os"}
    
    result = []
    for (element, comp_id), instances in groups.items():
        entry = dict(instances[0])  # first instance
        entry["count"] = len(instances)
        result.append(entry)

    result.sort(key=lambda m: (
        0 if m["element"] in PRIORITY_METALS else 1,
        m["count"],  # fewer instances = more interesting
        m["element"],
        m["comp_id"],
    ))

    return result

"""
Parser for .xyz coordinate files (DFT output format).

XYZ is the standard output from Gaussian, ORCA, ADF, etc.
Coordinates are already Cartesian — no cell parameters needed.

Format:
    N_atoms
    Comment line (often energy or method)
    Element  x  y  z
    Element  x  y  z
    ...
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from terrahedral.core import Atom, Bond, CoordinationSite, METAL_ELEMENTS
from terrahedral.transforms import distance


def parse_xyz(filepath: str | Path) -> list[dict]:
    """
    Parse a .xyz file.

    Returns list of dicts with element, x, y, z.
    Handles multi-frame trajectories (returns first frame).
    """
    lines = Path(filepath).read_text(errors="replace").strip().split("\n")
    if not lines:
        raise ValueError("Empty .xyz file")

    try:
        n_atoms = int(lines[0].strip())
    except ValueError:
        raise ValueError(f"First line of .xyz must be atom count, got: {lines[0]!r}")

    comment = lines[1].strip() if len(lines) > 1 else ""

    atoms = []
    for i in range(2, min(2 + n_atoms, len(lines))):
        parts = lines[i].split()
        if len(parts) < 4:
            continue
        try:
            atoms.append({
                "element": parts[0].capitalize(),
                "x": float(parts[1]),
                "y": float(parts[2]),
                "z": float(parts[3]),
            })
        except ValueError:
            continue

    return atoms


def load_xyz(
    filepath: str | Path,
    metal: Optional[str] = None,
    *,
    cutoff: float = 3.0,
) -> CoordinationSite:
    """
    Load a .xyz file and extract coordination around a metal center.

    Parameters
    ----------
    filepath : path to .xyz file
    metal : element symbol of the metal center (e.g. "Fe").
            If None, auto-detects the first metal in METAL_ELEMENTS.
    cutoff : distance cutoff for finding coordinated atoms (Å)

    Returns
    -------
    CoordinationSite centered on the metal.
    """
    atoms_raw = parse_xyz(filepath)

    if not atoms_raw:
        raise ValueError("No atoms found in .xyz file")

    # Find metal center
    metal_idx = None
    if metal:
        for i, a in enumerate(atoms_raw):
            if a["element"] == metal.capitalize():
                metal_idx = i
                break
        if metal_idx is None:
            available = [a["element"] for a in atoms_raw if a["element"] in METAL_ELEMENTS]
            raise ValueError(
                f"Metal '{metal}' not found. Available metals: {set(available)}"
            )
    else:
        # Auto-detect first metal
        for i, a in enumerate(atoms_raw):
            if a["element"] in METAL_ELEMENTS:
                metal_idx = i
                break
        if metal_idx is None:
            raise ValueError("No metal atom found in .xyz file. Specify metal= explicitly.")

    ma = atoms_raw[metal_idx]
    metal_xyz = (ma["x"], ma["y"], ma["z"])
    metal_atom = Atom(
        label=f"{ma['element']}1",
        element=ma["element"],
        x=0.0, y=0.0, z=0.0,
    )

    # Find coordinated atoms within cutoff
    ligand_atoms = []
    coord_bonds = []
    label_counter: dict[str, int] = {}

    for i, a in enumerate(atoms_raw):
        if i == metal_idx:
            continue
        # Skip H atoms (usually not direct ligands)
        if a["element"] == "H":
            continue
        xyz = (a["x"], a["y"], a["z"])
        d = distance(metal_xyz, xyz)
        if d <= cutoff:
            el = a["element"]
            label_counter[el] = label_counter.get(el, 0) + 1
            label = f"{el}{label_counter[el]}"

            ligand_atoms.append(Atom(
                label=label,
                element=el,
                x=a["x"] - ma["x"],
                y=a["y"] - ma["y"],
                z=a["z"] - ma["z"],
            ))
            coord_bonds.append(Bond(
                atom1=metal_atom.label,
                atom2=label,
                distance=d,
            ))

    return CoordinationSite(
        metal=metal_atom,
        ligands=ligand_atoms,
        bonds=coord_bonds,
        title=f"DFT — {ma['element']}",
        source=str(Path(filepath).stem),
    )

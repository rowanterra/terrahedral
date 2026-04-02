"""
Parser for standard PDB format files (ATOM/HETATM records).

Handles:
  - MD trajectory snapshots (.pdb from AMBER, GROMACS, etc.)
  - Deposited PDB structures downloaded directly
  - Single-model files (takes first MODEL if multiple)

Column layout (PDB fixed-width):
  1-6   Record type (ATOM/HETATM)
  7-11  Serial number
  13-16 Atom name
  17    Alternate location indicator
  18-20 Residue name
  22    Chain ID
  23-26 Residue sequence number
  31-38 X coordinate (Å)
  39-46 Y coordinate (Å)
  47-54 Z coordinate (Å)
  77-78 Element symbol
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from terrahedral.core import Atom, Bond, CoordinationSite, METAL_ELEMENTS
from terrahedral.transforms import distance


def parse_pdb(filepath: str | Path) -> list[dict]:
    """
    Parse a .pdb file and return list of atom dicts.

    Returns list of dicts with element, x, y, z, atom_name, resname, resseq, chain.
    Only reads the first MODEL if multiple models are present.
    """
    lines = Path(filepath).read_text(errors="replace").split("\n")

    atoms = []
    in_model = False

    for line in lines:
        rec = line[:6].strip()

        if rec == "MODEL":
            if in_model:
                break  # Stop at second model
            in_model = True
            continue

        if rec == "ENDMDL":
            break

        if rec not in ("ATOM", "HETATM"):
            continue

        try:
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21:22].strip()
            resseq = line[22:26].strip()
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])

            # Element symbol: columns 77-78 if available, else infer from atom name
            element = ""
            if len(line) > 77:
                raw_el = line[76:78].strip()
                # Reject placeholder values like "??" or empty
                if raw_el and raw_el != "??" and raw_el != "?":
                    element = raw_el

            if not element:
                # For standard residues (ATOM records), the element is
                # determined by the atom name. PDB convention:
                #   CA = C-alpha (Carbon), CB = C-beta, CD = C-delta,
                #   CE = C-epsilon, CG = C-gamma, CZ = C-zeta,
                #   NZ = N-zeta, NE = N-epsilon, ND = N-delta,
                #   OD = O-delta, OE = O-epsilon, OG = O-gamma,
                #   SD = S-delta, SG = S-gamma, HG = H-gamma, etc.
                # These are NOT Ca, Cb, Cd, Ce, Cg, Hg elements!
                #
                # Only HETATM records for standalone metal residues
                # (like "EU", "LA", "FE") use the atom name as element.
                PROTEIN_ATOMS = {
                    "C", "N", "O", "S", "H",
                    "CA", "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2",
                    "CE", "CE1", "CE2", "CE3", "CZ", "CZ2", "CZ3", "CH2",
                    "ND1", "ND2", "NE", "NE1", "NE2", "NZ", "NH1", "NH2",
                    "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH", "OXT",
                    "SD", "SG",
                }
                # Common H atom names (subset)
                is_hydrogen = atom_name.startswith("H") or atom_name in {
                    "HN", "HA", "HA2", "HA3", "HB", "HB1", "HB2", "HB3",
                    "HD1", "HD2", "HD3", "HD11", "HD12", "HD13",
                    "HD21", "HD22", "HD23",
                    "HE", "HE1", "HE2", "HE3", "HE21", "HE22",
                    "HG", "HG1", "HG2", "HG3", "HG11", "HG12", "HG13",
                    "HG21", "HG22", "HG23",
                    "HH", "HH11", "HH12", "HH21", "HH22",
                    "HZ", "HZ1", "HZ2", "HZ3",
                }

                if is_hydrogen:
                    element = "H"
                elif atom_name in PROTEIN_ATOMS:
                    element = atom_name[0]  # C, N, O, S
                elif rec == "ATOM":
                    # ATOM record: usually protein. But MD output
                    # sometimes puts metals as ATOM with element-like
                    # names (e.g. "LA" for lanthanum, "EU3" for europium).
                    name_clean = atom_name.strip().rstrip("0123456789+- ")
                    if (len(name_clean) <= 2 and len(name_clean) >= 1 and
                        name_clean.capitalize() in METAL_ELEMENTS and
                        resname.strip().rstrip("0123456789").upper() == name_clean.upper()):
                        # Atom name (stripped of digits) matches residue name
                        # and is a known metal — e.g. EU3/EU3, LA/LA
                        element = name_clean.capitalize()
                    else:
                        element = atom_name[0].upper()
                else:
                    # HETATM: could be a metal. Try stripping digits then
                    # matching 2-char or 1-char element symbols.
                    name_clean = atom_name.strip().rstrip("0123456789+- ")
                    if len(name_clean) >= 1 and len(name_clean) <= 2 and name_clean.capitalize() in METAL_ELEMENTS:
                        element = name_clean.capitalize()
                    elif atom_name.strip():
                        element = atom_name[0].upper()
                    else:
                        element = "?"

            atoms.append({
                "element": element.capitalize(),
                "x": x,
                "y": y,
                "z": z,
                "atom_name": atom_name,
                "resname": resname,
                "resseq": resseq,
                "chain": chain,
                "record": rec,
            })
        except (ValueError, IndexError):
            continue

    return atoms


def find_metals_pdb(filepath: str | Path) -> list[dict]:
    """Find all metal atoms in a PDB file."""
    atoms = parse_pdb(filepath)
    metals = []
    seen = set()

    for a in atoms:
        el = a["element"]
        if el in METAL_ELEMENTS:
            key = (el, a["resname"], a["resseq"], a["chain"])
            if key not in seen:
                seen.add(key)
                metals.append({
                    "element": el,
                    "label": f"{a['resname']}{a['resseq']}_{a['atom_name']}" if a["resname"] else f"{el}1",
                    "resname": a["resname"],
                    "resseq": a["resseq"],
                    "chain": a["chain"],
                    "x": a["x"],
                    "y": a["y"],
                    "z": a["z"],
                })

    return metals


def load_pdb(
    filepath: str | Path,
    metal: Optional[str] = None,
    *,
    cutoff: float = 3.0,
) -> CoordinationSite:
    """
    Load a .pdb file and extract coordination around a metal center.

    Parameters
    ----------
    filepath : path to .pdb file
    metal : element symbol or label of the metal center.
            If None, auto-detects the first metal in METAL_ELEMENTS.
    cutoff : distance cutoff for finding coordinated atoms (Å)

    Returns
    -------
    CoordinationSite centered on the metal.
    """
    atoms_raw = parse_pdb(filepath)

    if not atoms_raw:
        raise ValueError("No atoms found in .pdb file")

    # Lanthanides/actinides need wider cutoff
    LANTHANIDES = {"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"}
    ACTINIDES = {"Th","Pa","U","Np","Pu","Am"}

    # Find metal center
    metal_idx = None
    metal_el = None

    if metal:
        metal_cap = metal.capitalize()
        # Try element match first
        for i, a in enumerate(atoms_raw):
            if a["element"] == metal_cap:
                metal_idx = i
                metal_el = metal_cap
                break
        # Try label/atom_name match
        if metal_idx is None:
            for i, a in enumerate(atoms_raw):
                label = f"{a['resname']}{a['resseq']}_{a['atom_name']}"
                if label == metal or a["atom_name"] == metal.upper():
                    if a["element"] in METAL_ELEMENTS:
                        metal_idx = i
                        metal_el = a["element"]
                        break
    else:
        for i, a in enumerate(atoms_raw):
            if a["element"] in METAL_ELEMENTS:
                metal_idx = i
                metal_el = a["element"]
                break

    if metal_idx is None:
        available = sorted(set(a["element"] for a in atoms_raw if a["element"] in METAL_ELEMENTS))
        raise ValueError(
            f"Metal '{metal}' not found. Available metals: {available}"
        )

    # Widen cutoff for lanthanides/actinides
    if metal_el in LANTHANIDES or metal_el in ACTINIDES:
        cutoff = max(cutoff, 3.2)

    ma = atoms_raw[metal_idx]
    metal_xyz = (ma["x"], ma["y"], ma["z"])

    metal_label = ma["atom_name"] or f"{metal_el}1"
    metal_atom = Atom(
        label=metal_label,
        element=metal_el,
        x=0.0, y=0.0, z=0.0,
    )

    # Find coordinated atoms within cutoff (skip H)
    DONOR_ELEMENTS = {"N", "O", "S", "Se", "Cl", "Br", "I", "F"}

    ligand_atoms = []
    coord_bonds = []
    seen_labels: set[str] = set()

    candidates = []
    for i, a in enumerate(atoms_raw):
        if i == metal_idx:
            continue
        if a["element"] == "H":
            continue
        if a["element"] not in DONOR_ELEMENTS:
            continue

        xyz = (a["x"], a["y"], a["z"])
        d = distance(metal_xyz, xyz)
        if d <= cutoff and d > 0.5:
            candidates.append((d, i, a))

    candidates.sort(key=lambda x: x[0])

    # Limit CN based on metal type
    max_cn = 12 if (metal_el in LANTHANIDES or metal_el in ACTINIDES) else 7

    for d, idx, a in candidates[:max_cn]:
        # Build a unique label
        if a["resname"]:
            label = f"{a['resname']}{a['resseq']}_{a['atom_name']}"
        else:
            label = f"{a['element']}{len(ligand_atoms)+1}"

        # Deduplicate
        base_label = label
        suffix = 2
        while label in seen_labels:
            label = f"{base_label}_{suffix}"
            suffix += 1
        seen_labels.add(label)

        residue = f"{a['resname']}{a['resseq']}" if a["resname"] else ""

        ligand_atoms.append(Atom(
            label=label,
            element=a["element"],
            x=a["x"] - ma["x"],
            y=a["y"] - ma["y"],
            z=a["z"] - ma["z"],
            residue=residue,
            atom_name=a["atom_name"],
        ))
        coord_bonds.append(Bond(
            atom1=metal_atom.label,
            atom2=label,
            distance=round(d, 3),
        ))

    # Build title from filename
    stem = Path(filepath).stem
    title = f"{stem} — {metal_el}"

    return CoordinationSite(
        metal=metal_atom,
        ligands=ligand_atoms,
        bonds=coord_bonds,
        title=title,
        source=stem,
    )

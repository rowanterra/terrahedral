"""
High-level I/O: parse CIF files and return CoordinationSite objects.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional

from terrahedral.core import Atom, Bond, Angle, CoordinationSite, METAL_ELEMENTS
from terrahedral.transforms import frac_to_cart, distance, angle_at_center


# ── Constants shared across loaders ────────────────────────────────────────

# Donor elements that can coordinate to metals
DONOR_ELEMENTS = {"N", "O", "S", "Se", "Cl", "Br", "I", "F"}

# Backbone atom names that are typically NOT metal donors
BACKBONE_SKIP = {"C", "CA", "CB", "N", "H", "HA", "HB", "HN", "O"}

# Subset of BACKBONE_SKIP for metals that DO coordinate backbone carbonyl O
# (lanthanides, actinides, alkaline earths in EF-hand / calmodulin-like sites)
BACKBONE_SKIP_LANTHACTIN = {"C", "CA", "CB", "N", "H", "HA", "HB", "HN"}
# NOTE: backbone "O" (carbonyl) is intentionally ALLOWED for these metals

# Metals that commonly bind backbone carbonyl oxygens
BACKBONE_O_METALS = {
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",       # lanthanides
    "Th", "Pa", "U", "Np", "Pu", "Am",                # actinides
    "Ca", "Sr", "Ba",                                  # alkaline earths
}

# Lanthanide and actinide elements (need wider bond cutoffs)
LANTHANIDES = {
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd",
    "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
}
ACTINIDES = {"Th", "Pa", "U", "Np", "Pu", "Am"}


def load_cif(
    filepath: str | Path,
    metal: str,
    *,
    cutoff: float = 3.0,
) -> CoordinationSite:
    """
    Load a small-molecule CIF and extract coordination around *metal*.

    For mineral CIFs with symmetry operations, the asymmetric unit is
    expanded to a full unit cell (+ neighboring images) before finding
    coordination bonds.

    Parameters
    ----------
    filepath : path to .cif file
    metal : atom label of the metal center (e.g. "Pt9", "Cu1", "Fe")
    cutoff : max bond distance to consider (Angstrom), used as fallback
             if no _geom_bond records are present

    Returns
    -------
    CoordinationSite centered on the metal.
    """
    from terrahedral.parsers.cif_parser import parse, expand_symmetry

    cif = parse(filepath)
    cell = cif.cell

    # Symmetry expansion for mineral CIFs
    use_expansion = len(cif.symmetry_ops) > 1 and cell
    if use_expansion:
        expanded = expand_symmetry(cif)
    else:
        expanded = cif.atoms

    # Convert all atoms to Cartesian
    cart_lookup: dict[str, tuple[float, float, float]] = {}
    elem_lookup: dict[str, str] = {}

    for a in expanded:
        if a.get("element", "?") == "?":
            continue

        # Use pre-computed Cartesian if available (from expansion)
        if "cart_x" in a:
            xyz = (a["cart_x"], a["cart_y"], a["cart_z"])
        elif cell:
            xyz = frac_to_cart(
                a["frac_x"], a["frac_y"], a["frac_z"],
                cell["a"], cell["b"], cell["c"],
                cell["alpha"], cell["beta"], cell["gamma"],
            )
        else:
            xyz = (a["frac_x"], a["frac_y"], a["frac_z"])

        cart_lookup[a["label"]] = xyz
        elem_lookup[a["label"]] = a["element"]

    # Find the metal -- match by label first, then by element name
    if metal not in cart_lookup:
        element_matches = [
            lbl for lbl, el in elem_lookup.items() if el == metal
        ]
        if element_matches:
            metal = element_matches[0]
        else:
            available = [a["label"] for a in expanded if a["element"] in METAL_ELEMENTS]
            raise ValueError(
                f"Metal '{metal}' not found. Available metal labels: {available}"
            )

    mx, my, mz = cart_lookup[metal]

    # Find coordination bonds
    ligand_labels: set[str] = set()
    coord_bonds: list[Bond] = []

    if cif.bonds and not use_expansion:
        for b in cif.bonds:
            if b["atom1"] == metal:
                coord_bonds.append(Bond(metal, b["atom2"], b["distance"]))
                ligand_labels.add(b["atom2"])
            elif b["atom2"] == metal:
                coord_bonds.append(Bond(metal, b["atom1"], b["distance"]))
                ligand_labels.add(b["atom1"])
    else:
        # Distance-based neighbor detection
        metal_el = elem_lookup[metal]
        min_d = 1.5  # Minimum realistic bond distance

        all_neighbors = []
        for label, xyz in cart_lookup.items():
            if label == metal:
                continue
            d = distance(cart_lookup[metal], xyz)
            if d <= cutoff and d > min_d:
                lig_el = elem_lookup.get(label, "?")
                all_neighbors.append((label, d, lig_el))

        if use_expansion:
            # Mineral mode: build coordination shell intelligently
            # Sort by distance, prefer non-same-element neighbors
            all_neighbors.sort(key=lambda x: x[1])

            # First pass: find the natural gap in bond distances
            # (coordination shell vs next shell)
            if all_neighbors:
                dists = [d for _, d, _ in all_neighbors]
                # Find biggest gap in first 12 neighbors
                gap_idx = len(dists)
                for i in range(1, min(12, len(dists))):
                    if dists[i] - dists[i-1] > 0.4:  # 0.4 Å gap
                        gap_idx = i
                        break

                # Take neighbors up to the gap, or all within cutoff
                for label, d, lig_el in all_neighbors[:max(gap_idx, 2)]:
                    coord_bonds.append(Bond(metal, label, round(d, 3)))
                    ligand_labels.add(label)
        else:
            for label, d, lig_el in all_neighbors:
                coord_bonds.append(Bond(metal, label, round(d, 3)))
                ligand_labels.add(label)

    # Find angles at the metal
    coord_angles: list[Angle] = []
    for a in cif.angles:
        if a["center"] == metal:
            coord_angles.append(Angle(a["atom1"], a["center"], a["atom2"], a["angle"]))

    # Build Atom objects, centered on metal
    metal_atom = Atom(
        label=metal,
        element=elem_lookup[metal],
        x=0.0, y=0.0, z=0.0,
    )

    ligand_atoms = []
    for lbl in sorted(ligand_labels):
        if lbl in cart_lookup:
            x, y, z = cart_lookup[lbl]
            ligand_atoms.append(Atom(
                label=lbl,
                element=elem_lookup.get(lbl, "?"),
                x=x - mx, y=y - my, z=z - mz,
            ))

    return CoordinationSite(
        metal=metal_atom,
        ligands=ligand_atoms,
        bonds=coord_bonds,
        angles=coord_angles,
        title=f"{cif.formula} — {elem_lookup.get(metal, metal)}",
        source=str(Path(filepath).stem),
    )


def load_mmcif(
    filepath: str | Path,
    metal_asym: str,
    metal_seq: str = ".",
    *,
    site_label: Optional[str] = None,
    metal_element: Optional[str] = None,
    shell_mode: str = "auto",
) -> CoordinationSite:
    """
    Load a macromolecular mmCIF and extract coordination for a specific metal.

    Parameters
    ----------
    filepath : path to PDB mmCIF
    metal_asym : asym_id (chain) of the metal HETATM
    metal_seq : seq_id of the metal ("." for single-atom HETATM)
    site_label : optional display label (e.g. "T1Cu")
    metal_element : optional element symbol to disambiguate metals in same
                    residue (e.g. "Mn" vs "Ca" in OEX cluster)
    shell_mode : coordination shell detection strategy. One of:

        ``"auto"`` (default)
            Adaptive behaviour.  When PDB ``_struct_conn`` metalc records
            list fewer than 4 bonds for a transition-metal site, the gap
            filter is *skipped* so that supplement-found ligands are kept
            (addresses incomplete depositor annotations at lower resolution).
            When metalc records already list ≥ 6 bonds, the supplement
            search is *capped* at +1 extra atom so that crystallisation
            additives and second-shell contacts are less likely to inflate
            the CN.  Otherwise the standard gap filter is applied.

        ``"strict"``
            Always apply the gap-based shell filter with the original
            thresholds (MIN_GAP = 0.35 Å, ≥ 3 retained, post-gap > 2.5 Å)
            and allow unlimited supplement additions.  Best for high-
            throughput scans where over-counting is the bigger risk.

        ``"permissive"``
            Never apply the gap filter.  All atoms from metalc records
            *and* the 3.5 Å supplement search are retained.  Best for
            curated analyses where every candidate ligand should be
            inspected manually.

    Returns
    -------
    CoordinationSite centered on the metal.
    """
    if shell_mode not in ("auto", "strict", "permissive"):
        raise ValueError(
            f"shell_mode must be 'auto', 'strict', or 'permissive', "
            f"got {shell_mode!r}"
        )
    from terrahedral.parsers.mmcif_parser import parse, METAL_RESIDUES

    mmcif = parse(filepath)

    # Find the metal atom coordinates
    _metal_xyz = None
    _metal_element = None
    _metal_comp = None
    for a in mmcif.atoms:
        if a["asym_id"] == metal_asym:
            seq_match = (metal_seq == "." or a["seq_id"] == metal_seq)
            is_metal = (a["comp_id"] in METAL_RESIDUES or
                        a.get("element", "") in METAL_ELEMENTS)
            # If caller specified an element, require exact match
            if metal_element and a.get("element", "") != metal_element:
                continue
            if seq_match and is_metal:
                _metal_xyz = (a["x"], a["y"], a["z"])
                _metal_element = a["element"]
                _metal_comp = a["comp_id"]
                break

    # Assign to the names used by the rest of the function
    metal_xyz = _metal_xyz
    metal_element_found = _metal_element
    metal_comp = _metal_comp

    if metal_xyz is None:
        from terrahedral.parsers.mmcif_parser import find_metals
        metals = find_metals(mmcif)
        labels = [f"{m['asym_id']}/{m['seq_id']} ({m['comp_id']})" for m in metals]
        raise ValueError(
            f"Metal at asym={metal_asym} seq={metal_seq} not found. "
            f"Available metals: {labels}"
        )

    mx, my, mz = metal_xyz

    # Find metalc bonds to this metal
    coord_bonds: list[Bond] = []
    ligand_infos: list[dict] = []

    for b in mmcif.metalc_bonds:
        # Determine which partner is the metal and which is the ligand
        # Match by comp_id (from METAL_RESIDUES) or by the specific metal comp we found
        ptnr2_is_metal = (b["ptnr2_comp"] in METAL_RESIDUES or
                          b["ptnr2_comp"] == metal_comp)
        ptnr1_is_metal = (b["ptnr1_comp"] in METAL_RESIDUES or
                          b["ptnr1_comp"] == metal_comp)

        if ptnr2_is_metal and b["ptnr2_asym"] == metal_asym:
            ligand_infos.append({
                "asym": b["ptnr1_asym"],
                "comp": b["ptnr1_comp"],
                "seq": b["ptnr1_seq"],
                "atom": b["ptnr1_atom"],
                "alt": b.get("ptnr1_alt", ""),
                "distance": b["distance"],
                "auth_seq": b.get("auth1_seq", ""),
            })
        elif ptnr1_is_metal and b["ptnr1_asym"] == metal_asym:
            ligand_infos.append({
                "asym": b["ptnr2_asym"],
                "comp": b["ptnr2_comp"],
                "seq": b["ptnr2_seq"],
                "atom": b["ptnr2_atom"],
                "alt": b.get("ptnr2_alt", ""),
                "distance": b["distance"],
                "auth_seq": b.get("auth2_seq", ""),
            })

    # Resolve ligand atom coordinates
    metal_label = site_label or f"{metal_comp}{metal_seq}"
    metal_atom = Atom(
        label=metal_label,
        element=metal_element_found,
        x=0.0, y=0.0, z=0.0,
    )

    ligand_atoms: list[Atom] = []
    seen_labels: set[str] = set()
    symop_warnings: list[str] = []

    for lig in ligand_infos:
        if lig["distance"] is None:
            continue

        # Find matching atom in coordinate records
        for a in mmcif.atoms:
            if (
                a["comp_id"] == lig["comp"]
                and a["atom_name"] == lig["atom"]
                and a["asym_id"] == lig["asym"]
            ):
                # Match on seq_id if available and specific
                if lig["seq"] and lig["seq"] != "." and a["seq_id"] != lig["seq"]:
                    continue

                # For HOH/water and other HETATM with seq_id=".", disambiguate
                # by auth_seq_id — critical when multiple waters exist
                if lig["auth_seq"] and a.get("auth_seq_id", ""):
                    if a["auth_seq_id"] != lig["auth_seq"]:
                        continue

                # Handle alt conformers: take first match or matching alt
                if lig["alt"] and a.get("alt_id", "") and a["alt_id"] != lig["alt"]:
                    continue

                auth_seq = a.get("auth_seq_id", lig["seq"])
                label = f"{lig['comp']}{auth_seq}_{lig['atom']}"

                if label in seen_labels:
                    continue

                # Distance validation: reject symmetry-mate atoms whose
                # deposited coordinates don't match the CIF bond distance.
                # Symmetry operators like 3_675 mean the real ligand is in
                # a different unit cell copy — we only have identity coords.
                cart_d = distance((a["x"], a["y"], a["z"]), metal_xyz)
                if lig["distance"] and abs(cart_d - lig["distance"]) > 0.5:
                    # Wrong symmetry copy — flag and skip
                    symop_warnings.append(
                        f"{lig['comp']}{auth_seq}.{lig['atom']} "
                        f"(cart={cart_d:.1f} vs cif={lig['distance']:.1f} Å, "
                        f"likely symmetry mate)"
                    )
                    continue

                seen_labels.add(label)

                residue = f"{lig['comp']}{auth_seq}"

                ligand_atoms.append(Atom(
                    label=label,
                    element=a["element"],
                    x=a["x"] - mx,
                    y=a["y"] - my,
                    z=a["z"] - mz,
                    residue=residue,
                    atom_name=lig["atom"],
                ))

                coord_bonds.append(Bond(
                    atom1=metal_label,
                    atom2=label,
                    distance=lig["distance"],
                ))
                break

    if symop_warnings:
        for w in symop_warnings:
            warnings.warn(f"Skipped symmetry-mate ligand: {w}", stacklevel=2)

    # ── Supplementary neighbor search ──
    # _struct_conn metalc records are often incomplete — PDB frequently
    # omits weak/long axial bonds (e.g. Met-S(thioether) at ~2.9 Å in
    # Type 1 blue copper sites like rusticyanin/plastocyanin).
    # When metalc records found *some* ligands, scan for additional
    # close donor atoms that were missed.
    # PDB _struct_conn metalc records frequently undercount for
    # lanthanides (often listing 7-8 of 9 bonds in EF-hand sites).
    if ligand_atoms:
        # Use a wider cutoff for lanthanides/actinides where CN 8-10
        # is common and metalc annotations are often incomplete.
        if metal_element_found in LANTHANIDES or metal_element_found in ACTINIDES:
            supplement_cutoff = 3.5
        else:
            supplement_cutoff = 3.5

        # Labels already found via metalc
        existing_labels = seen_labels.copy()

        for a in mmcif.atoms:
            if a["element"] not in DONOR_ELEMENTS:
                continue

            ax, ay, az = a["x"], a["y"], a["z"]
            d = distance((ax, ay, az), metal_xyz)
            if d > supplement_cutoff or d < 0.5:
                continue

            # Skip atoms already picked up from metalc records
            auth_seq = a.get("auth_seq_id", a.get("seq_id", ""))
            comp = a.get("comp_id", "")
            atom_name = a.get("atom_name", "")
            label = f"{comp}{auth_seq}_{atom_name}"
            if label in existing_labels:
                continue

            record = a.get("record", "ATOM")

            # Only consider sidechain donor atoms (SD, SG, OD1, OD2,
            # OE1, OE2, NE2, ND1, NZ, OH, OG, OG1, etc.) and HETATM
            # For lanthanides/actinides/Ca, backbone carbonyl O is a
            # legitimate donor (EF-hand / calmodulin-like sites).
            skip_set = (BACKBONE_SKIP_LANTHACTIN
                        if metal_element_found in BACKBONE_O_METALS
                        else BACKBONE_SKIP)
            if record == "ATOM" and atom_name in skip_set:
                continue
            if a["element"] == "H":
                continue

            # Must be within the same chain or a HETATM
            # (avoid picking up donors from neighboring chains)
            if record == "ATOM" and a["asym_id"] != metal_asym:
                # Allow if it's the polymer chain that the metal is
                # coordinated to (check existing ligands for chain)
                if a["asym_id"] not in {li_info["asym"] for li_info in ligand_infos}:
                    continue

            existing_labels.add(label)
            residue = f"{comp}{auth_seq}"

            ligand_atoms.append(Atom(
                label=label,
                element=a["element"],
                x=ax - mx,
                y=ay - my,
                z=az - mz,
                residue=residue,
                atom_name=atom_name,
            ))

            coord_bonds.append(Bond(
                atom1=metal_label,
                atom2=label,
                distance=round(d, 3),
                is_long=d > 2.6,
            ))

    # ── Distance-based fallback ──
    # If no ligands were found from _struct_conn metalc records,
    # search for nearby donor atoms within a cutoff.
    # Common for lanthanides, actinides, and entries with missing annotations.
    if not ligand_atoms:
        # Cutoff depends on metal type — lanthanides/actinides have longer bonds
        if metal_element_found in LANTHANIDES or metal_element_found in ACTINIDES:
            cutoff = 3.2
        elif metal_element_found in ("Ca","Sr","Ba","K","Na","Rb","Cs"):
            cutoff = 3.2
        else:
            cutoff = 2.8

        candidates = []
        for a in mmcif.atoms:
            if a["element"] not in DONOR_ELEMENTS:
                continue
            if a["element"] in METAL_ELEMENTS:
                continue

            ax, ay, az = a["x"], a["y"], a["z"]
            d = distance((ax, ay, az), metal_xyz)
            if d > cutoff or d < 0.5:  # too far or too close (same atom)
                continue

            atom_name = a.get("atom_name", "")
            record = a.get("record", "ATOM")
            comp_id = a.get("comp_id", "")

            # Skip backbone atoms from protein residues (ATOM records)
            # but allow sidechain donors and HETATM donors
            # For lanthanides/actinides/Ca, backbone carbonyl O is allowed.
            # Note: backbone N *can* be a donor in some cases — allow it
            # only if the atom is from a HETATM or is very close
            skip_set_fb = (BACKBONE_SKIP_LANTHACTIN
                           if metal_element_found in BACKBONE_O_METALS
                           else BACKBONE_SKIP)
            if record == "ATOM" and atom_name in skip_set_fb:
                # Exception: allow backbone N only if very close (< 2.5 Å)
                if atom_name == "N" and d < 2.5:
                    pass
                else:
                    continue

            # Skip hydrogen
            if a["element"] == "H":
                continue

            candidates.append({
                "atom": a,
                "distance": round(d, 3),
                "atom_name": atom_name,
                "comp_id": comp_id,
            })

        # Sort by distance, take closest neighbors up to reasonable CN
        candidates.sort(key=lambda c: c["distance"])
        # Limit: typically CN ≤ 10 for lanthanides, ≤ 6 for transition metals
        max_cn = 10 if (metal_element_found in LANTHANIDES or metal_element_found in ACTINIDES) else 7

        seen_fallback: set[str] = set()
        for cand in candidates[:max_cn]:
            a = cand["atom"]
            auth_seq = a.get("auth_seq_id", a.get("seq_id", ""))
            comp = a.get("comp_id", "")
            atom_name = cand["atom_name"]

            label = f"{comp}{auth_seq}_{atom_name}"

            # Deduplicate
            if label in seen_fallback:
                suffix = 2
                while f"{label}_{suffix}" in seen_fallback:
                    suffix += 1
                label = f"{label}_{suffix}"
            seen_fallback.add(label)

            residue = f"{comp}{auth_seq}"

            ligand_atoms.append(Atom(
                label=label,
                element=a["element"],
                x=a["x"] - mx,
                y=a["y"] - my,
                z=a["z"] - mz,
                residue=residue,
                atom_name=atom_name,
            ))

            coord_bonds.append(Bond(
                atom1=metal_label,
                atom2=label,
                distance=cand["distance"],
            ))

    # ── Gap-based coordination shell filter ──
    # After collecting all candidate ligands (from metalc + supplement +
    # fallback), trim second-sphere atoms by finding a significant distance
    # gap.  This prevents inflated CN from PDB metalc records that list
    # second-shell interactions (e.g. 2.7–2.9 Å for Fe³⁺ whose real
    # first-shell Fe-O bonds are 1.9–2.2 Å).
    #
    # The filter is controlled by shell_mode:
    #   "strict"     — always apply (original behaviour)
    #   "permissive" — never apply (keep everything from metalc + supplement)
    #   "auto"       — skip when metalc records appear incomplete (< 4 bonds
    #                  for a transition metal), apply otherwise

    # Count how many ligands came from the authoritative metalc records
    # (i.e. those that were present *before* the supplement search).
    n_metalc = len(ligand_infos)  # from _struct_conn parsing above

    apply_gap_filter = True
    if shell_mode == "permissive":
        apply_gap_filter = False
    elif shell_mode == "auto":
        # Heuristic: if _struct_conn listed fewer than 5 bonds for a
        # transition metal, the annotations are likely incomplete and the
        # supplement-found atoms should be trusted rather than trimmed.
        # Threshold 5 chosen because common TM geometries (octahedral,
        # trigonal bipyramidal) have CN ≥ 5; sparse metalc + gap filter
        # often incorrectly trims real ligands that sit slightly further.
        is_transition = (metal_element_found not in LANTHANIDES
                         and metal_element_found not in ACTINIDES
                         and metal_element_found not in ("Ca", "Sr", "Ba",
                                                         "K", "Na", "Rb",
                                                         "Cs", "Mg"))
        if is_transition and n_metalc < 5:
            apply_gap_filter = False

    if apply_gap_filter and len(ligand_atoms) > 1:
        # Pair ligands with their bond distances
        lig_dist = []
        dist_map = {b.atom2: b.distance for b in coord_bonds if b.distance}
        for la in ligand_atoms:
            d = dist_map.get(la.label)
            if d is None:
                d = distance((la.x, la.y, la.z), (0, 0, 0))
            lig_dist.append((la, d))
        lig_dist.sort(key=lambda x: x[1])

        # Find the largest gap between consecutive distances
        distances = [d for _, d in lig_dist]
        best_gap = 0.0
        best_gap_idx = len(distances) - 1  # default: keep all

        for i in range(len(distances) - 1):
            gap = distances[i + 1] - distances[i]
            if gap > best_gap:
                best_gap = gap
                best_gap_idx = i

        # Only trim if:
        # - the gap is significant (> 0.35 Å)
        # - it leaves at least 3 ligands (avoid cutting to trivial CN)
        # - the shell beyond the gap is clearly second-sphere (> 2.5 Å)
        MIN_GAP = 0.35
        if (best_gap > MIN_GAP
                and best_gap_idx + 1 >= 3
                and distances[best_gap_idx + 1] > 2.5):
            keep_labels = {la.label for la, _ in lig_dist[:best_gap_idx + 1]}
            ligand_atoms = [la for la in ligand_atoms if la.label in keep_labels]
            coord_bonds = [b for b in coord_bonds
                           if b.atom2 in keep_labels or b.atom1 in keep_labels]

    # ── Auto-mode supplement cap ──
    # When metalc annotations are already rich (≥ 6 bonds) and we're in
    # auto mode, cap the total CN to n_metalc to prevent crystallisation
    # additives or second-shell contacts from inflating the count.
    # Rationale: if the depositor annotated ≥ 6 bonds the coordination
    # shell is likely complete; supplement-found atoms beyond that are
    # usually PEG oxygens, buffer molecules, or second-shell contacts.
    if shell_mode == "auto" and n_metalc >= 6 and len(ligand_atoms) > n_metalc:
        # Keep the n_metalc closest atoms
        lig_dist = []
        dist_map = {b.atom2: b.distance for b in coord_bonds if b.distance}
        for la in ligand_atoms:
            d = dist_map.get(la.label)
            if d is None:
                d = distance((la.x, la.y, la.z), (0, 0, 0))
            lig_dist.append((la, d))
        lig_dist.sort(key=lambda x: x[1])
        keep_labels = {la.label for la, _ in lig_dist[:n_metalc]}
        ligand_atoms = [la for la in ligand_atoms if la.label in keep_labels]
        coord_bonds = [b for b in coord_bonds
                       if b.atom2 in keep_labels or b.atom1 in keep_labels]

    site = CoordinationSite(
        metal=metal_atom,
        ligands=ligand_atoms,
        bonds=coord_bonds,
        title=site_label or mmcif.title,
        source=str(Path(filepath).stem).upper(),
        resolution=mmcif.resolution,
    )
    site._symop_warnings = list(symop_warnings)
    return site

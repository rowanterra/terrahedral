"""
Core data structures for coordination environments.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence


# ── Atom colors by element ──────────────────────────────────────────────────

ELEMENT_STYLE = {
    # Transition metals — first row
    "Ti": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "V":  {"fill": "#E8593C", "stroke": "#A93220", "r": 22, "txt": "#fff"},
    "Cr": {"fill": "#BA7517", "stroke": "#854F0B", "r": 22, "txt": "#fff"},
    "Mn": {"fill": "#D4537E", "stroke": "#993556", "r": 22, "txt": "#fff"},
    "Fe": {"fill": "#D85A30", "stroke": "#993C1D", "r": 22, "txt": "#fff"},
    "Co": {"fill": "#1D9E75", "stroke": "#0F6E56", "r": 22, "txt": "#fff"},
    "Ni": {"fill": "#97C459", "stroke": "#3B6D11", "r": 22, "txt": "#fff"},
    "Cu": {"fill": "#378ADD", "stroke": "#185FA5", "r": 22, "txt": "#fff"},
    "Zn": {"fill": "#7F77DD", "stroke": "#534AB7", "r": 22, "txt": "#fff"},
    # Transition metals — second row
    "Zr": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "Nb": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "Mo": {"fill": "#7F77DD", "stroke": "#534AB7", "r": 22, "txt": "#fff"},
    "Tc": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "Ru": {"fill": "#1D9E75", "stroke": "#0F6E56", "r": 22, "txt": "#fff"},
    "Rh": {"fill": "#D85A30", "stroke": "#993C1D", "r": 22, "txt": "#fff"},
    "Pd": {"fill": "#888780", "stroke": "#5F5E5A", "r": 22, "txt": "#fff"},
    "Ag": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 22, "txt": "#2C2C2A"},
    "Cd": {"fill": "#EF9F27", "stroke": "#854F0B", "r": 22, "txt": "#fff"},
    # Transition metals — third row
    "Hf": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "Ta": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "W":  {"fill": "#888780", "stroke": "#5F5E5A", "r": 22, "txt": "#fff"},
    "Re": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    "Os": {"fill": "#378ADD", "stroke": "#185FA5", "r": 22, "txt": "#fff"},
    "Ir": {"fill": "#1D9E75", "stroke": "#0F6E56", "r": 22, "txt": "#fff"},
    "Pt": {"fill": "#888780", "stroke": "#5F5E5A", "r": 24, "txt": "#fff"},
    "Au": {"fill": "#EFB827", "stroke": "#A97F0B", "r": 22, "txt": "#412402"},
    "Hg": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 22, "txt": "#2C2C2A"},
    # Lanthanides (shared teal-violet palette)
    "La": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Ce": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Pr": {"fill": "#7BAEBE", "stroke": "#4A8A99", "r": 24, "txt": "#fff"},
    "Nd": {"fill": "#8BB8C8", "stroke": "#5A9AA9", "r": 24, "txt": "#fff"},
    "Pm": {"fill": "#8BB8C8", "stroke": "#5A9AA9", "r": 24, "txt": "#fff"},
    "Sm": {"fill": "#8BB8C8", "stroke": "#5A9AA9", "r": 24, "txt": "#fff"},
    "Eu": {"fill": "#9BC2D2", "stroke": "#6AAAB9", "r": 24, "txt": "#fff"},
    "Gd": {"fill": "#9BC2D2", "stroke": "#6AAAB9", "r": 24, "txt": "#fff"},
    "Tb": {"fill": "#ABCCD8", "stroke": "#7ABAC9", "r": 24, "txt": "#2C2C2A"},
    "Dy": {"fill": "#ABCCD8", "stroke": "#7ABAC9", "r": 24, "txt": "#2C2C2A"},
    "Ho": {"fill": "#BBD6DE", "stroke": "#8ACAD9", "r": 24, "txt": "#2C2C2A"},
    "Er": {"fill": "#BBD6DE", "stroke": "#8ACAD9", "r": 24, "txt": "#2C2C2A"},
    "Tm": {"fill": "#CBE0E8", "stroke": "#9ADAE9", "r": 24, "txt": "#2C2C2A"},
    "Yb": {"fill": "#CBE0E8", "stroke": "#9ADAE9", "r": 24, "txt": "#2C2C2A"},
    "Lu": {"fill": "#DBE8EE", "stroke": "#AAE0E9", "r": 24, "txt": "#2C2C2A"},
    "Y":  {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 22, "txt": "#fff"},
    "Sc": {"fill": "#7A8B99", "stroke": "#4A5B69", "r": 22, "txt": "#fff"},
    # Alkaline / alkaline earth
    "Li": {"fill": "#9575CD", "stroke": "#6A3FAA", "r": 14, "txt": "#fff"},
    "Be": {"fill": "#66BB6A", "stroke": "#338836", "r": 14, "txt": "#fff"},
    "Mg": {"fill": "#66BB6A", "stroke": "#338836", "r": 18, "txt": "#fff"},
    "Ca": {"fill": "#66BB6A", "stroke": "#338836", "r": 20, "txt": "#fff"},
    "Na": {"fill": "#9575CD", "stroke": "#6A3FAA", "r": 18, "txt": "#fff"},
    "K":  {"fill": "#9575CD", "stroke": "#6A3FAA", "r": 20, "txt": "#fff"},
    "Rb": {"fill": "#9575CD", "stroke": "#6A3FAA", "r": 22, "txt": "#fff"},
    "Cs": {"fill": "#9575CD", "stroke": "#6A3FAA", "r": 24, "txt": "#fff"},
    "Sr": {"fill": "#66BB6A", "stroke": "#338836", "r": 22, "txt": "#fff"},
    "Ba": {"fill": "#66BB6A", "stroke": "#338836", "r": 24, "txt": "#fff"},
    # Post-transition / metalloids
    "Al": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 18, "txt": "#2C2C2A"},
    "Ga": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 18, "txt": "#2C2C2A"},
    "Ge": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 18, "txt": "#2C2C2A"},
    "In": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 20, "txt": "#2C2C2A"},
    "Sn": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 20, "txt": "#2C2C2A"},
    "Tl": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 22, "txt": "#2C2C2A"},
    "Pb": {"fill": "#7F77DD", "stroke": "#534AB7", "r": 22, "txt": "#fff"},
    "Bi": {"fill": "#D4537E", "stroke": "#993556", "r": 22, "txt": "#fff"},
    "Sb": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 20, "txt": "#2C2C2A"},
    "Si": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 16, "txt": "#2C2C2A"},
    "As": {"fill": "#B4B2A9", "stroke": "#7A7875", "r": 16, "txt": "#2C2C2A"},
    "Te": {"fill": "#EF9F27", "stroke": "#854F0B", "r": 16, "txt": "#412402"},
    # Actinides
    "Th": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Pa": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "U":  {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Np": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Pu": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    "Am": {"fill": "#6BA3BE", "stroke": "#3A7A99", "r": 24, "txt": "#fff"},
    # Donor atoms
    "N":  {"fill": "#5DCAA5", "stroke": "#0F6E56", "r": 14, "txt": "#fff"},
    "O":  {"fill": "#F09595", "stroke": "#C03030", "r": 14, "txt": "#fff"},
    "S":  {"fill": "#FAC775", "stroke": "#854F0B", "r": 14, "txt": "#412402"},
    "P":  {"fill": "#ED93B1", "stroke": "#993556", "r": 16, "txt": "#4B1528"},
    "B":  {"fill": "#97C459", "stroke": "#3B6D11", "r": 14, "txt": "#173404"},
    "C":  {"fill": "#B4B2A9", "stroke": "#5F5E5A", "r": 12, "txt": "#2C2C2A"},
    "Cl": {"fill": "#97C459", "stroke": "#3B6D11", "r": 14, "txt": "#173404"},
    "Br": {"fill": "#D85A30", "stroke": "#993C1D", "r": 14, "txt": "#fff"},
    "Se": {"fill": "#EF9F27", "stroke": "#854F0B", "r": 14, "txt": "#412402"},
    "I":  {"fill": "#7F77DD", "stroke": "#534AB7", "r": 14, "txt": "#fff"},
    "F":  {"fill": "#90EE90", "stroke": "#3B8E3B", "r": 12, "txt": "#173404"},
}

_DEFAULT_STYLE = {"fill": "#B4B2A9", "stroke": "#5F5E5A", "r": 14, "txt": "#2C2C2A"}

# Comprehensive set of metal elements for detection in CIF files.
METAL_ELEMENTS = {
    # Alkali / alkaline earth
    "Li", "Na", "K", "Rb", "Cs", "Be", "Mg", "Ca", "Sr", "Ba",
    # First-row transition
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    # Second-row transition
    "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
    # Third-row transition
    "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    # Post-transition / metalloids (appear as metal-like HETATM in PDB)
    "Al", "Ga", "Ge", "In", "Sn", "Tl", "Pb", "Bi", "Sb",
    "As", "Se", "Te",
    # Lanthanides (all 15)
    "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu",
    # Actinides
    "Th", "Pa", "U", "Np", "Pu", "Am",
}


# ── Data classes ────────────────────────────────────────────────────────────

@dataclass
class Atom:
    """Single atom with Cartesian coordinates and display metadata."""

    label: str
    element: str
    x: float
    y: float
    z: float
    residue: str = ""
    atom_name: str = ""

    @property
    def xyz(self) -> tuple[float, float, float]:
        return (self.x, self.y, self.z)

    @property
    def style(self) -> dict:
        return ELEMENT_STYLE.get(self.element, _DEFAULT_STYLE)

    def distance_to(self, other: Atom) -> float:
        return math.sqrt(
            (self.x - other.x) ** 2
            + (self.y - other.y) ** 2
            + (self.z - other.z) ** 2
        )

    def centered_on(self, origin: Atom) -> Atom:
        """Return copy translated so *origin* is at (0, 0, 0)."""
        return Atom(
            label=self.label,
            element=self.element,
            x=self.x - origin.x,
            y=self.y - origin.y,
            z=self.z - origin.z,
            residue=self.residue,
            atom_name=self.atom_name,
        )


@dataclass
class Bond:
    """Metal-ligand coordination bond."""

    atom1: str          # label of first atom (usually the metal)
    atom2: str          # label of second atom (ligand donor)
    distance: float     # Angstrom
    bond_type: str = "coordination"
    is_long: bool = False   # e.g. weak axial bond

    @property
    def label(self) -> str:
        return f"{self.atom1}-{self.atom2}"


@dataclass
class Angle:
    """Ligand-metal-ligand angle."""

    atom1: str
    center: str
    atom2: str
    value: float        # degrees

    @property
    def label(self) -> str:
        return f"{self.atom1}-{self.center}-{self.atom2}"


# ── Coordination site ───────────────────────────────────────────────────────

@dataclass
class CoordinationSite:
    """
    A metal coordination environment: metal center + first-shell ligands
    + bond distances + angles.

    This is the primary object you work with. Create from a CIF file::

        site = CoordinationSite.from_cif("structure.cif", metal="Pt9")
        site = CoordinationSite.from_mmcif("6GTL.cif", metal_asym="B", metal_seq=".")

    Or build manually::

        site = CoordinationSite(
            metal=Atom("Cu", "Cu", 0, 0, 0),
            ligands=[...],
            bonds=[...],
        )
    """

    metal: Atom
    ligands: list[Atom] = field(default_factory=list)
    bonds: list[Bond] = field(default_factory=list)
    angles: list[Angle] = field(default_factory=list)
    title: str = ""
    source: str = ""          # e.g. "PDB 6GTL" or "COD 4070506"
    geometry: str = ""        # e.g. "distorted tetrahedral"
    resolution: float = 0.0   # Angstrom, from diffraction data
    space_group: str = ""
    formula: str = ""
    _symop_warnings: list[str] = field(default_factory=list, repr=False)

    # ── Properties ──

    @property
    def coordination_number(self) -> int:
        return len(self.ligands)

    @property
    def donor_set(self) -> str:
        """Sorted donor atom formula, e.g. 'N2S2' or 'His3NO2'."""
        elements = sorted(lig.element for lig in self.ligands)
        result = []
        i = 0
        while i < len(elements):
            el = elements[i]
            count = 1
            while i + count < len(elements) and elements[i + count] == el:
                count += 1
            result.append(f"{el}{count}" if count > 1 else el)
            i += count
        return "".join(result)

    @property
    def bond_dict(self) -> dict[str, float]:
        """Map ligand label -> bond distance."""
        return {b.atom2: b.distance for b in self.bonds}

    # ── Constructors ──

    @classmethod
    def from_cif(
        cls,
        filepath: str | Path,
        metal: str,
        *,
        cutoff: float = 3.0,
    ) -> CoordinationSite:
        """
        Build from a small-molecule CIF (COD/CCDC format).

        Parameters
        ----------
        filepath : path to .cif file
        metal : atom label of the metal center (e.g. "Pt9", "Cu1")
        cutoff : max bond distance to consider (Angstrom), used as fallback
                 if no _geom_bond records are present
        """
        from terrahedral.io import load_cif
        return load_cif(filepath, metal, cutoff=cutoff)

    @classmethod
    def from_mmcif(
        cls,
        filepath: str | Path,
        metal_asym: str,
        metal_seq: str = ".",
        *,
        site_label: Optional[str] = None,
        metal_element: Optional[str] = None,
        shell_mode: str = "auto",
    ) -> CoordinationSite:
        """
        Build from a macromolecular mmCIF (PDB format).

        Parameters
        ----------
        filepath : path to .cif file from PDB
        metal_asym : chain / asym_id containing the metal HETATM
        metal_seq : seq_id (residue number) of the metal, or "." for HETATM
        site_label : optional label like "T1Cu" for the title
        metal_element : optional element symbol to disambiguate
        shell_mode : coordination shell detection strategy
        """
        from terrahedral.io import load_mmcif
        return load_mmcif(filepath, metal_asym, metal_seq,
                          site_label=site_label, metal_element=metal_element,
                          shell_mode=shell_mode)

    @classmethod
    def from_pdb(
        cls,
        pdb_id: str,
        metal_asym: str,
        metal_seq: str = ".",
        *,
        site_label: Optional[str] = None,
        cache_dir: Optional[str | Path] = None,
        force: bool = False,
    ) -> CoordinationSite:
        """
        Fetch an mmCIF from RCSB PDB and build a CoordinationSite.

        Parameters
        ----------
        pdb_id : 4-character PDB identifier (e.g. "6GTL")
        metal_asym : chain / asym_id containing the metal HETATM
        metal_seq : seq_id of the metal, or "." for single-atom HETATM
        site_label : optional label like "T1Cu"
        cache_dir : where to store downloaded files
        force : re-download even if cached
        """
        from terrahedral.fetch import fetch_pdb
        path = fetch_pdb(pdb_id, cache_dir=cache_dir, force=force)
        return cls.from_mmcif(path, metal_asym, metal_seq, site_label=site_label)

    @classmethod
    def from_cod(
        cls,
        cod_id: str,
        metal: str,
        *,
        cutoff: float = 3.0,
        cache_dir: Optional[str | Path] = None,
        force: bool = False,
    ) -> CoordinationSite:
        """
        Fetch a CIF from the Crystallography Open Database and build a CoordinationSite.

        Parameters
        ----------
        cod_id : numeric COD identifier (e.g. "4070506")
        metal : atom label of the metal center
        cutoff : max bond distance fallback (Angstrom)
        cache_dir : where to store downloaded files
        force : re-download even if cached
        """
        from terrahedral.fetch import fetch_cod
        path = fetch_cod(cod_id, cache_dir=cache_dir, force=force)
        return cls.from_cif(path, metal, cutoff=cutoff)

    @classmethod
    def from_xyz(
        cls,
        filepath: str | Path,
        metal: Optional[str] = None,
        *,
        cutoff: float = 3.0,
    ) -> CoordinationSite:
        """
        Build from a DFT .xyz coordinate file.

        Parameters
        ----------
        filepath : path to .xyz file
        metal : element symbol to center on (auto-detects if None)
        cutoff : max bond distance (Angstrom)
        """
        from terrahedral.parsers.xyz_parser import load_xyz
        return load_xyz(filepath, metal=metal, cutoff=cutoff)

    # ── Analysis shortcuts ──

    def classify_geometry(self) -> str:
        """Assign geometry label (e.g. 'square pyramidal') via tau analysis."""
        from terrahedral.analysis.geometry import classify
        self.geometry = classify(self)
        return self.geometry

    def tau5(self) -> Optional[float]:
        """Addison tau-5 parameter for 5-coordinate sites."""
        from terrahedral.analysis.geometry import tau5
        return tau5(self)

    def tau4(self) -> Optional[float]:
        """Yang tau-4 parameter for 4-coordinate sites."""
        from terrahedral.analysis.geometry import tau4
        return tau4(self)

    def ideal_geometry(self) -> Optional[str]:
        """Return the nearest ideal geometry, or None if already ideal."""
        from terrahedral.analysis import ideal_geometry
        return ideal_geometry(self)

    def deviation_summary(self) -> str:
        """Describe how geometry deviates from ideal (tau values, classification)."""
        from terrahedral.analysis import deviation_summary
        return deviation_summary(self)

    # ── Display ──

    def __repr__(self) -> str:
        geo = f" {self.geometry}" if self.geometry else ""
        return (
            f"CoordinationSite({self.metal.element}, "
            f"CN={self.coordination_number}, "
            f"donors={self.donor_set}{geo})"
        )

    def summary(self) -> str:
        """Human-readable text summary."""
        lines = [
            f"{'─' * 50}",
            f"  {self.title or self.metal.label}" + (f"  ({self.source})" if self.source else ""),
            f"  Metal: {self.metal.element} at ({self.metal.x:.3f}, {self.metal.y:.3f}, {self.metal.z:.3f})",
            f"  CN = {self.coordination_number}  Donors: {self.donor_set}",
        ]
        if self.geometry:
            lines.append(f"  Geometry: {self.geometry}")
        lines.append(f"  {'─' * 46}")
        lines.append(f"  {'Ligand':<20s} {'Element':<8s} {'Distance (A)':<14s}")
        lines.append(f"  {'─' * 46}")
        for lig in self.ligands:
            d = self.bond_dict.get(lig.label, lig.distance_to(self.metal))
            res = f" ({lig.residue})" if lig.residue else ""
            lines.append(f"  {lig.label + res:<20s} {lig.element:<8s} {d:<14.3f}")
        if self.angles:
            lines.append(f"  {'─' * 46}")
            lines.append(f"  Key angles:")
            for a in sorted(self.angles, key=lambda a: -a.value)[:6]:
                lines.append(f"    {a.label}: {a.value:.1f} deg")
        lines.append(f"{'─' * 50}")
        return "\n".join(lines)

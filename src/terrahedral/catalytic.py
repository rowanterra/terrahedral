"""
CatalyticFrame — represents a catalytic intermediate or state
constructed from spectroscopic observables, DFT coordinates,
or crystallographic data.

Designed for the Guo-lab workflow:
    spectroscopy (Mössbauer, EPR, NRVS) → DFT → coordinates → geometry

Usage:

    # From spectroscopic observables (no coordinates needed)
    frame = CatalyticFrame.manual(
        metal="Fe", oxidation_state=4,
        ligands=["N", "N", "O", "O", "O"],
        bond_lengths=[2.10, 2.05, 1.65, 2.00, 2.12],
        mossbauer={"delta": 0.10, "delta_eq": 1.05},
        epr={"spin": 1, "g": [2.27, 2.21, 1.98]},
        label="FeIV=O intermediate",
    )
    frame.site.render_svg(style="orbital", oxidation_state=4)

    # From DFT .xyz file
    frame = CatalyticFrame.from_xyz("dft_optimized.xyz", metal="Fe",
                                     oxidation_state=4, spin=1)

    # Build a reaction pathway from multiple frames
    pathway = CatalyticPathway([frame1, frame2, frame3])
    pathway.render_svg("pathway.svg")
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Sequence

from terrahedral.core import Atom, Bond, CoordinationSite


# ─────────────────────────────────────────────────────────────────────
# Spectroscopic parameter lookup tables
# ─────────────────────────────────────────────────────────────────────

# Mössbauer isomer shift (δ, mm/s) ranges for Fe oxidation states + spin
# at 4.2K vs α-Fe.  Ranges are approximate — overlap exists.
MOSSBAUER_DELTA_RANGES: dict[tuple[str, int, str], tuple[float, float]] = {
    # (element, oxidation_state, spin) → (delta_min, delta_max)
    # ── Iron ──
    ("Fe", 2, "high"):  (0.80, 1.40),   # HS Fe²⁺
    ("Fe", 2, "low"):   (0.20, 0.50),   # LS Fe²⁺
    ("Fe", 3, "high"):  (0.30, 0.55),   # HS Fe³⁺
    ("Fe", 3, "low"):   (0.10, 0.30),   # LS Fe³⁺
    ("Fe", 4, "high"):  (0.00, 0.20),   # HS Fe⁴⁺ (rare)
    ("Fe", 4, "low"):   (-0.10, 0.15),  # Fe(IV)=O ferryl
    ("Fe", 5, "low"):   (-0.10, 0.05),  # Fe(V) — very rare
    ("Fe", 6, "low"):   (-0.10, 0.00),  # Ferrate(VI)
}

# Quadrupole splitting (ΔEQ) rough ranges (mm/s)
MOSSBAUER_DEQ_RANGES: dict[tuple[str, int, str], tuple[float, float]] = {
    ("Fe", 2, "high"):  (1.50, 3.50),
    ("Fe", 2, "low"):   (0.20, 1.00),
    ("Fe", 3, "high"):  (0.30, 1.00),
    ("Fe", 3, "low"):   (1.50, 3.00),
    ("Fe", 4, "high"):  (0.50, 1.50),
    ("Fe", 4, "low"):   (0.80, 2.00),   # Ferryl typically 0.9-1.7
}

# NRVS / vibrational frequency → bond length correlations
# frequency_cm (cm⁻¹) → approximate bond_length (Å)
NRVS_CORRELATIONS: dict[str, list[tuple[float, float, float, float]]] = {
    # bond_type: [(freq_min, freq_max, length_min, length_max), ...]
    "Fe=O":  [(790, 870, 1.62, 1.70)],    # Ferryl oxo stretch
    "Fe-O":  [(400, 600, 1.85, 2.20)],    # Fe–OH or Fe–O(carboxylate)
    "Fe-OH": [(550, 650, 1.80, 1.95)],    # Fe–OH hydroxide
    "Fe-N":  [(350, 500, 2.00, 2.25)],    # Fe–N(His/imidazole)
    "Fe-S":  [(300, 400, 2.25, 2.45)],    # Fe–S(Cys)
    "Mn=O":  [(600, 700, 1.65, 1.80)],    # Mn(IV/V)=O
    "Cu-N":  [(350, 450, 1.95, 2.10)],    # Cu–N(His)
    "Cu-S":  [(300, 380, 2.10, 2.30)],    # Cu–S(Cys/Met)
}


def interpret_mossbauer(
    delta: float,
    delta_eq: float,
    element: str = "Fe",
) -> list[dict]:
    """
    Interpret Mössbauer parameters.

    Returns list of possible assignments ranked by fit,
    each with oxidation_state, spin, confidence, notes.
    """
    candidates = []
    for (el, ox, spin), (d_min, d_max) in MOSSBAUER_DELTA_RANGES.items():
        if el != element:
            continue
        if d_min <= delta <= d_max:
            # Check ΔEQ if available
            deq_range = MOSSBAUER_DEQ_RANGES.get((el, ox, spin))
            deq_ok = True
            deq_note = ""
            if deq_range:
                deq_ok = deq_range[0] - 0.3 <= delta_eq <= deq_range[1] + 0.3
                if not deq_ok:
                    deq_note = f" (ΔEQ={delta_eq} outside typical {deq_range[0]}-{deq_range[1]})"

            # Score: how centered is delta in the range?
            center = (d_min + d_max) / 2
            width = d_max - d_min
            confidence = max(0, 1.0 - abs(delta - center) / (width / 2))
            if not deq_ok:
                confidence *= 0.5

            spin_label = "HS" if spin == "high" else "LS"
            candidates.append({
                "oxidation_state": ox,
                "spin": spin,
                "spin_label": spin_label,
                "confidence": round(confidence, 2),
                "notes": f"{element}({ox}+) {spin_label}, δ={delta} ΔEQ={delta_eq}{deq_note}",
            })

    return sorted(candidates, key=lambda c: -c["confidence"])


def freq_to_bond_length(
    frequency_cm: float,
    bond_type: str,
) -> Optional[float]:
    """
    Estimate bond length from vibrational frequency (NRVS/IR/Raman).

    Parameters
    ----------
    frequency_cm : vibrational frequency in cm⁻¹
    bond_type : e.g. "Fe=O", "Fe-N", "Fe-S"

    Returns
    -------
    Estimated bond length in Å, or None if out of range.
    """
    corrs = NRVS_CORRELATIONS.get(bond_type, [])
    for freq_min, freq_max, len_min, len_max in corrs:
        if freq_min <= frequency_cm <= freq_max:
            # Linear interpolation (higher freq → shorter bond)
            frac = (frequency_cm - freq_min) / (freq_max - freq_min)
            return len_max - frac * (len_max - len_min)
    return None


# ─────────────────────────────────────────────────────────────────────
# CatalyticFrame
# ─────────────────────────────────────────────────────────────────────

@dataclass
class CatalyticFrame:
    """
    A single catalytic intermediate or state.

    Wraps a CoordinationSite with spectroscopic metadata:
    oxidation state, spin state, spectroscopic observables,
    and an interpretive label.
    """

    site: CoordinationSite
    oxidation_state: int = 2
    spin: Optional[int] = None          # total spin quantum number S
    spin_label: str = ""                # "HS" / "LS" / "S=1" etc.
    label: str = ""                     # e.g. "FeIV=O intermediate"
    step: int = 0                       # position in reaction pathway

    # Spectroscopic observables (if known)
    mossbauer: Optional[dict] = None    # {"delta": ..., "delta_eq": ...}
    epr: Optional[dict] = None          # {"spin": ..., "g": [...]}
    nrvs: Optional[dict] = None         # {"Fe=O": 831, "Fe-N": 420, ...}
    dft_energy: Optional[float] = None  # in kcal/mol or Hartree

    # Interpretation
    assignments: list = field(default_factory=list)

    @classmethod
    def manual(
        cls,
        metal: str,
        oxidation_state: int,
        ligands: list[str],
        bond_lengths: list[float],
        *,
        spin: Optional[int] = None,
        mossbauer: Optional[dict] = None,
        epr: Optional[dict] = None,
        nrvs: Optional[dict] = None,
        dft_energy: Optional[float] = None,
        label: str = "",
        step: int = 0,
        residues: Optional[list[str]] = None,
    ) -> CatalyticFrame:
        """
        Build a CatalyticFrame from spectroscopic observables.

        No coordinates or CIF file needed — constructs a synthetic
        CoordinationSite from the donor atoms and bond lengths.

        Parameters
        ----------
        metal : element symbol (e.g. "Fe")
        oxidation_state : metal oxidation state (e.g. 4)
        ligands : list of donor element symbols ["N", "N", "O", "O", "O"]
        bond_lengths : corresponding M–L distances in Å
        spin : total spin S (e.g. 1 for S=1 triplet)
        mossbauer : {"delta": float, "delta_eq": float}
        epr : {"spin": int, "g": [gx, gy, gz]}
        nrvs : {"Fe=O": 831, "Fe-N": 420} — frequencies in cm⁻¹
        label : human name like "FeIV=O intermediate"
        step : position in a reaction pathway (0, 1, 2, ...)
        residues : optional residue labels ["His214", "His268", ...]
        """
        if len(ligands) != len(bond_lengths):
            raise ValueError(
                f"ligands ({len(ligands)}) and bond_lengths ({len(bond_lengths)}) "
                f"must have the same length"
            )

        metal_atom = Atom(label=f"{metal}1", element=metal, x=0.0, y=0.0, z=0.0)

        # Distribute ligands evenly in 3D
        n = len(ligands)
        ligand_atoms = []
        bonds = []
        label_counter: dict[str, int] = {}

        for i, (el, d) in enumerate(zip(ligands, bond_lengths)):
            label_counter[el] = label_counter.get(el, 0) + 1
            lbl = f"{el}{label_counter[el]}"

            # Distribute on a sphere — uses golden spiral for even spacing
            theta = math.acos(1 - 2 * (i + 0.5) / n)
            phi = math.pi * (1 + 5 ** 0.5) * i
            x = d * math.sin(theta) * math.cos(phi)
            y = d * math.sin(theta) * math.sin(phi)
            z = d * math.cos(theta)

            res = residues[i] if residues and i < len(residues) else ""
            ligand_atoms.append(Atom(
                label=lbl, element=el, x=x, y=y, z=z, residue=res,
            ))
            bonds.append(Bond(
                atom1=metal_atom.label, atom2=lbl, distance=d,
                bond_type="oxo" if d < 1.75 else "coordination",
                is_long=d > 2.4,
            ))

        title = label or f"{metal}({oxidation_state}+) CN={n}"
        site = CoordinationSite(
            metal=metal_atom,
            ligands=ligand_atoms,
            bonds=bonds,
            title=title,
        )
        site.classify_geometry()

        # Determine spin label
        spin_label = ""
        if spin is not None:
            spin_label = f"S={spin}"
        elif epr and "spin" in epr:
            spin = epr["spin"]
            spin_label = f"S={spin}"

        # Interpret Mössbauer if provided
        assignments = []
        if mossbauer:
            assignments = interpret_mossbauer(
                mossbauer.get("delta", 0),
                mossbauer.get("delta_eq", 0),
                element=metal,
            )

        return cls(
            site=site,
            oxidation_state=oxidation_state,
            spin=spin,
            spin_label=spin_label,
            label=label,
            step=step,
            mossbauer=mossbauer,
            epr=epr,
            nrvs=nrvs,
            dft_energy=dft_energy,
            assignments=assignments,
        )

    @classmethod
    def from_xyz(
        cls,
        filepath: str | Path,
        metal: Optional[str] = None,
        *,
        oxidation_state: int = 2,
        spin: Optional[int] = None,
        cutoff: float = 3.0,
        label: str = "",
        step: int = 0,
        mossbauer: Optional[dict] = None,
        epr: Optional[dict] = None,
        nrvs: Optional[dict] = None,
        dft_energy: Optional[float] = None,
    ) -> CatalyticFrame:
        """
        Build a CatalyticFrame from a DFT .xyz file.

        Parameters
        ----------
        filepath : path to .xyz file
        metal : element to center on (auto-detects if None)
        oxidation_state : metal oxidation state
        spin : total spin S
        cutoff : bond distance cutoff in Å
        label : human name for this intermediate
        """
        from terrahedral.parsers.xyz_parser import load_xyz

        site = load_xyz(filepath, metal=metal, cutoff=cutoff)
        site.classify_geometry()

        spin_label = f"S={spin}" if spin is not None else ""

        assignments = []
        if mossbauer:
            assignments = interpret_mossbauer(
                mossbauer.get("delta", 0),
                mossbauer.get("delta_eq", 0),
                element=site.metal.element,
            )

        return cls(
            site=site,
            oxidation_state=oxidation_state,
            spin=spin,
            spin_label=spin_label,
            label=label or site.title,
            step=step,
            mossbauer=mossbauer,
            epr=epr,
            nrvs=nrvs,
            dft_energy=dft_energy,
            assignments=assignments,
        )

    @classmethod
    def from_site(
        cls,
        site: CoordinationSite,
        *,
        oxidation_state: int = 2,
        spin: Optional[int] = None,
        label: str = "",
        step: int = 0,
        **kwargs,
    ) -> CatalyticFrame:
        """Wrap an existing CoordinationSite as a CatalyticFrame."""
        return cls(
            site=site,
            oxidation_state=oxidation_state,
            spin=spin,
            spin_label=f"S={spin}" if spin is not None else "",
            label=label or site.title,
            step=step,
            **kwargs,
        )

    def summary(self) -> str:
        """Human-readable summary of this catalytic state."""
        lines = [f"{'━' * 50}"]
        lines.append(f"  {self.label or 'Catalytic Frame'}")
        lines.append(f"  {self.site.metal.element}({self.oxidation_state}+) "
                      f"CN={self.site.coordination_number} "
                      f"{self.spin_label} "
                      f"Geometry: {self.site.geometry}")

        if self.mossbauer:
            d = self.mossbauer.get("delta", "?")
            deq = self.mossbauer.get("delta_eq", "?")
            lines.append(f"  Mössbauer: δ = {d} mm/s, ΔEQ = {deq} mm/s")

        if self.epr:
            g = self.epr.get("g", [])
            lines.append(f"  EPR: S={self.epr.get('spin', '?')}, "
                         f"g = [{', '.join(f'{x:.2f}' for x in g)}]")

        if self.nrvs:
            for bond_type, freq in self.nrvs.items():
                est = freq_to_bond_length(freq, bond_type)
                est_str = f" → ~{est:.2f} Å" if est else ""
                lines.append(f"  NRVS: {bond_type} = {freq} cm⁻¹{est_str}")

        if self.dft_energy is not None:
            lines.append(f"  DFT energy: {self.dft_energy:.2f} kcal/mol")

        if self.assignments:
            lines.append(f"  Mössbauer assignment: {self.assignments[0]['notes']} "
                         f"(confidence: {self.assignments[0]['confidence']:.0%})")

        lines.append(f"  {'─' * 46}")
        for lig in self.site.ligands:
            d = self.site.bond_dict.get(lig.label, lig.distance_to(self.site.metal))
            res = f" ({lig.residue})" if lig.residue else ""
            btype = ""
            for b in self.site.bonds:
                if b.atom2 == lig.label:
                    if b.bond_type == "oxo":
                        btype = " [oxo]"
                    elif b.is_long:
                        btype = " [long]"
            lines.append(f"  {lig.label + res:<20s} {lig.element:<4s} {d:.3f} Å{btype}")

        lines.append(f"{'━' * 50}")
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────
# CatalyticPathway — a sequence of frames
# ─────────────────────────────────────────────────────────────────────

@dataclass
class CatalyticPathway:
    """
    A reaction pathway: an ordered sequence of CatalyticFrames
    representing intermediates in a catalytic cycle.
    """

    frames: list[CatalyticFrame] = field(default_factory=list)
    title: str = ""

    def add(self, frame: CatalyticFrame) -> None:
        """Add a frame. Auto-assigns step number."""
        frame.step = len(self.frames)
        self.frames.append(frame)

    def summary(self) -> str:
        """Text summary of the entire pathway."""
        lines = [f"{'═' * 60}"]
        lines.append(f"  Catalytic Pathway: {self.title or '(untitled)'}")
        lines.append(f"  {len(self.frames)} intermediates")
        lines.append(f"{'═' * 60}")

        for i, f in enumerate(self.frames):
            arrow = " → " if i > 0 else "   "
            geo = f.site.geometry or "?"
            energy = f"  ΔE={f.dft_energy:.1f}" if f.dft_energy is not None else ""
            lines.append(
                f"{arrow}[{i}] {f.label or f.site.metal.element}"
                f"({f.oxidation_state}+) {f.spin_label} "
                f"CN={f.site.coordination_number} {geo}{energy}"
            )

        lines.append(f"{'═' * 60}")
        return "\n".join(lines)

    def render_svg(
        self,
        output: Optional[str | Path] = None,
        width: int = 900,
    ) -> str:
        """
        Render pathway as a horizontal SVG showing each intermediate
        with arrows between them.
        """
        from terrahedral.core import ELEMENT_STYLE, _DEFAULT_STYLE

        n = len(self.frames)
        if n == 0:
            return "<svg></svg>"

        frame_w = max(140, (width - 40) // n)
        height = 320
        lines = [
            f'<svg width="100%" viewBox="0 0 {width} {height}" '
            f'xmlns="http://www.w3.org/2000/svg">',
            '<style>',
            '.pf-label{font-family:monospace;font-size:11px;'
            'fill:var(--color-text-primary,#222);text-anchor:middle}',
            '.pf-sub{font-family:monospace;font-size:9px;'
            'fill:var(--color-text-secondary,#666);text-anchor:middle}',
            '.pf-energy{font-family:monospace;font-size:9px;'
            'fill:var(--color-text-info,#185FA5);text-anchor:middle;font-weight:500}',
            '.pf-arrow{stroke:var(--color-text-tertiary,#999);stroke-width:1.5;'
            'fill:none;marker-end:url(#pfarr)}',
            '</style>',
            '<defs><marker id="pfarr" viewBox="0 0 8 8" refX="7" refY="4" '
            'markerWidth="5" markerHeight="5" orient="auto">'
            '<path d="M1 1L7 4L1 7" fill="none" '
            'stroke="var(--color-text-tertiary,#999)" stroke-width="1.2"/>'
            '</marker></defs>',
        ]

        # Title
        if self.title:
            lines.append(f'<text x="{width // 2}" y="24" class="pf-label" '
                         f'font-size="14" font-weight="600">{self.title}</text>')

        base_y = 160
        for i, frame in enumerate(self.frames):
            cx = 30 + i * frame_w + frame_w // 2

            # Arrow from previous
            if i > 0:
                ax1 = cx - frame_w + 30
                ax2 = cx - 30
                lines.append(f'<line x1="{ax1}" y1="{base_y}" '
                             f'x2="{ax2}" y2="{base_y}" class="pf-arrow"/>')

            # Metal circle
            ms = ELEMENT_STYLE.get(frame.site.metal.element, _DEFAULT_STYLE)
            r = 22
            lines.append(f'<circle cx="{cx}" cy="{base_y}" r="{r}" '
                         f'fill="{ms["fill"]}" stroke="{ms["stroke"]}" stroke-width="2"/>')
            lines.append(f'<text x="{cx}" y="{base_y + 4}" class="pf-label" '
                         f'fill="{ms["txt"]}" font-weight="600">'
                         f'{frame.site.metal.element}</text>')

            # Mini ligands
            nl = len(frame.site.ligands)
            for j, lig in enumerate(frame.site.ligands):
                angle = 2 * math.pi * j / max(nl, 1) - math.pi / 2
                lr = 38
                lx = cx + lr * math.cos(angle)
                ly = base_y + lr * math.sin(angle)
                ls = ELEMENT_STYLE.get(lig.element, _DEFAULT_STYLE)
                # Bond line
                is_oxo = any(b.atom2 == lig.label and b.bond_type == "oxo"
                             for b in frame.site.bonds)
                sw = "2.5" if is_oxo else "1.5"
                lines.append(f'<line x1="{cx}" y1="{base_y}" x2="{lx}" y2="{ly}" '
                             f'stroke="#999" stroke-width="{sw}"/>')
                lines.append(f'<circle cx="{lx}" cy="{ly}" r="8" '
                             f'fill="{ls["fill"]}" stroke="{ls["stroke"]}" stroke-width="1"/>')
                lines.append(f'<text x="{lx}" y="{ly + 3}" class="pf-sub" '
                             f'fill="{ls["txt"]}" font-size="7">{lig.element}</text>')

            # Labels below
            label = frame.label or f"Step {i}"
            lines.append(f'<text x="{cx}" y="{base_y + 60}" class="pf-label">'
                         f'{label}</text>')

            ox_spin = f"{frame.site.metal.element}({frame.oxidation_state}+)"
            if frame.spin_label:
                ox_spin += f" {frame.spin_label}"
            lines.append(f'<text x="{cx}" y="{base_y + 74}" class="pf-sub">'
                         f'{ox_spin}</text>')

            geo = frame.site.geometry or ""
            lines.append(f'<text x="{cx}" y="{base_y + 87}" class="pf-sub">'
                         f'{geo}</text>')

            if frame.dft_energy is not None:
                lines.append(f'<text x="{cx}" y="{base_y - 35}" class="pf-energy">'
                             f'ΔE = {frame.dft_energy:.1f}</text>')

        lines.append("</svg>")
        svg = "\n".join(lines)

        if output:
            Path(output).write_text(svg)
        return svg

"""
Electronic configurations, d- and f-orbital splitting data for all metals.

Maps metal + oxidation state → electron count, common geometries,
spectroscopic properties.  Covers d-block (transition metals),
f-block (lanthanides, actinides), and s-block (alkali/alkaline earth).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional


@dataclass
class ElectronicConfig:
    """Electronic configuration for a metal ion."""
    element: str
    oxidation_state: int
    d_count: int
    spin_state: str = "low"     # "low" or "high"
    color_origin: str = ""      # e.g. "d-d", "LMCT", "MLCT"
    notes: str = ""
    f_count: int = 0            # for lanthanides / actinides
    block: str = "d"            # "s", "d", or "f"


# ─────────────────────────────────────────────────────────────────────
# Electronic configurations — every metal at common oxidation states
# ─────────────────────────────────────────────────────────────────────

CONFIGS: dict[tuple[str, int], ElectronicConfig] = {
    # ══════════════════════════════════════════════════════════════════
    # s-block — alkali / alkaline earth
    # ══════════════════════════════════════════════════════════════════
    ("Li", 1): ElectronicConfig("Li", 1, 0, "low", "none",
                                "noble gas core, small hard Lewis acid", block="s"),
    ("Be", 2): ElectronicConfig("Be", 2, 0, "low", "none",
                                "no d electrons, toxic, tetrahedral preference", block="s"),
    ("Na", 1): ElectronicConfig("Na", 1, 0, "low", "none",
                                "no d electrons, weak electrostatic coordination", block="s"),
    ("K", 1):  ElectronicConfig("K",  1, 0, "low", "none",
                                "no d electrons, large ionic radius", block="s"),
    ("Rb", 1): ElectronicConfig("Rb", 1, 0, "low", "none",
                                "very large radius, weak coordination", block="s"),
    ("Cs", 1): ElectronicConfig("Cs", 1, 0, "low", "none",
                                "largest common monocation", block="s"),
    ("Mg", 2): ElectronicConfig("Mg", 2, 0, "low", "none",
                                "d⁰: no d electrons, electrostatic binding only", block="s"),
    ("Ca", 2): ElectronicConfig("Ca", 2, 0, "low", "none",
                                "d⁰: no d electrons, structural / signaling", block="s"),
    ("Sr", 2): ElectronicConfig("Sr", 2, 0, "low", "none",
                                "Ca mimic, larger radius", block="s"),
    ("Ba", 2): ElectronicConfig("Ba", 2, 0, "low", "none",
                                "very large radius, high CN preferred", block="s"),

    # ══════════════════════════════════════════════════════════════════
    # d-block — first-row transition metals
    # ══════════════════════════════════════════════════════════════════
    # Scandium
    ("Sc", 3): ElectronicConfig("Sc", 3, 0, "low", "none",
                                "d⁰: smallest trivalent group 3 ion"),
    # Titanium
    ("Ti", 2): ElectronicConfig("Ti", 2, 2, "high", "d-d",
                                "d²: rare, strong reductant"),
    ("Ti", 3): ElectronicConfig("Ti", 3, 1, "low", "d-d",
                                "d¹: purple in aqueous, EPR active"),
    ("Ti", 4): ElectronicConfig("Ti", 4, 0, "low", "LMCT",
                                "d⁰: colorless, TiO₂ white pigment"),
    # Vanadium
    ("V", 2):  ElectronicConfig("V",  2, 3, "high", "d-d",
                                "d³: violet, strong reductant"),
    ("V", 3):  ElectronicConfig("V",  3, 2, "low", "d-d",
                                "d²: green in octahedral"),
    ("V", 4):  ElectronicConfig("V",  4, 1, "low", "d-d",
                                "d¹: one unpaired electron, EPR active"),
    ("V", 5):  ElectronicConfig("V",  5, 0, "low", "LMCT",
                                "d⁰: no d-d, color from O→V LMCT"),
    # Chromium
    ("Cr", 2): ElectronicConfig("Cr", 2, 4, "high", "d-d",
                                "d⁴: Jahn-Teller, sky blue, strong reductant"),
    ("Cr", 3): ElectronicConfig("Cr", 3, 3, "high", "d-d",
                                "d³: half-filled t₂g, kinetically inert"),
    ("Cr", 6): ElectronicConfig("Cr", 6, 0, "low", "LMCT",
                                "d⁰: chromate yellow from LMCT"),
    # Manganese
    ("Mn", 2): ElectronicConfig("Mn", 2, 5, "high", "d-d (weak)",
                                "d⁵ HS: all d-d transitions spin-forbidden"),
    ("Mn", 3): ElectronicConfig("Mn", 3, 4, "high", "d-d / LMCT",
                                "d⁴: Jahn-Teller active in octahedral"),
    ("Mn", 4): ElectronicConfig("Mn", 4, 3, "high", "d-d",
                                "d³: half-filled t₂g in octahedral"),
    ("Mn", 7): ElectronicConfig("Mn", 7, 0, "low", "LMCT",
                                "d⁰: permanganate purple from LMCT"),
    # Iron
    ("Fe", 2): ElectronicConfig("Fe", 2, 6, "variable", "d-d",
                                "d⁶: LS = diamagnetic, HS = 4 unpaired"),
    ("Fe", 3): ElectronicConfig("Fe", 3, 5, "high", "d-d / LMCT",
                                "d⁵: spin-forbidden d-d, weak color"),
    ("Fe", 4): ElectronicConfig("Fe", 4, 4, "high", "d-d",
                                "d⁴: found in ferryl (FeIV=O) intermediates"),
    # Cobalt
    ("Co", 2): ElectronicConfig("Co", 2, 7, "variable", "d-d",
                                "d⁷: HS tetrahedral intensely colored"),
    ("Co", 3): ElectronicConfig("Co", 3, 6, "low", "d-d",
                                "d⁶ LS: diamagnetic, kinetically inert"),
    # Nickel
    ("Ni", 2): ElectronicConfig("Ni", 2, 8, "variable", "d-d",
                                "d⁸: square planar = diamagnetic"),
    ("Ni", 3): ElectronicConfig("Ni", 3, 7, "low", "d-d",
                                "d⁷: found in Ni-SOD, Jahn-Teller"),
    # Copper
    ("Cu", 1): ElectronicConfig("Cu", 1, 10, "low", "none (d¹⁰)",
                                "d¹⁰: no d-d transitions, colorless"),
    ("Cu", 2): ElectronicConfig("Cu", 2, 9, "low", "d-d / LMCT",
                                "d⁹: Jahn-Teller active in octahedral"),
    # Zinc
    ("Zn", 2): ElectronicConfig("Zn", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰: spectroscopically silent"),

    # ══════════════════════════════════════════════════════════════════
    # d-block — second-row transition metals
    # ══════════════════════════════════════════════════════════════════
    ("Y", 3):  ElectronicConfig("Y",  3, 0, "low", "none",
                                "d⁰: pseudo-lanthanide, diamagnetic"),
    ("Zr", 4): ElectronicConfig("Zr", 4, 0, "low", "LMCT",
                                "d⁰: hard Lewis acid, oxophilic"),
    ("Nb", 5): ElectronicConfig("Nb", 5, 0, "low", "LMCT",
                                "d⁰: isoelectronic with ZrIV"),
    ("Nb", 3): ElectronicConfig("Nb", 3, 2, "low", "d-d",
                                "d²: strong-field, paired"),
    ("Mo", 4): ElectronicConfig("Mo", 4, 2, "low", "d-d",
                                "d²: diamagnetic in strong-field"),
    ("Mo", 5): ElectronicConfig("Mo", 5, 1, "low", "d-d",
                                "d¹: EPR active, blue/green"),
    ("Mo", 6): ElectronicConfig("Mo", 6, 0, "low", "LMCT",
                                "d⁰: common in oxo-transfer enzymes"),
    ("Tc", 7): ElectronicConfig("Tc", 7, 0, "low", "LMCT",
                                "d⁰: pertechnetate, medical imaging"),
    ("Tc", 4): ElectronicConfig("Tc", 4, 3, "low", "d-d",
                                "d³: common in ⁹⁹ᵐTc complexes"),
    ("Ru", 2): ElectronicConfig("Ru", 2, 6, "low", "d-d / MLCT",
                                "d⁶ LS: diamagnetic, strong-field"),
    ("Ru", 3): ElectronicConfig("Ru", 3, 5, "low", "d-d / LMCT",
                                "d⁵: one unpaired (LS), paramagnetic"),
    ("Rh", 1): ElectronicConfig("Rh", 1, 8, "low", "d-d",
                                "d⁸: square planar, Wilkinson's catalyst"),
    ("Rh", 3): ElectronicConfig("Rh", 3, 6, "low", "d-d",
                                "d⁶ LS: octahedral, inert"),
    ("Pd", 2): ElectronicConfig("Pd", 2, 8, "low", "d-d",
                                "d⁸: square planar, cross-coupling"),
    ("Pd", 4): ElectronicConfig("Pd", 4, 6, "low", "d-d",
                                "d⁶: octahedral, strong oxidant"),
    ("Ag", 1): ElectronicConfig("Ag", 1, 10, "low", "none (d¹⁰)",
                                "d¹⁰: colorless, linear preference"),
    ("Ag", 2): ElectronicConfig("Ag", 2, 9, "low", "d-d",
                                "d⁹: rare, strong oxidant, paramagnetic"),
    ("Cd", 2): ElectronicConfig("Cd", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰: spectroscopically silent, Zn mimic"),

    # ══════════════════════════════════════════════════════════════════
    # d-block — third-row transition metals
    # ══════════════════════════════════════════════════════════════════
    ("Hf", 4): ElectronicConfig("Hf", 4, 0, "low", "LMCT",
                                "d⁰: isoelectronic with ZrIV, hard acid"),
    ("Ta", 5): ElectronicConfig("Ta", 5, 0, "low", "LMCT",
                                "d⁰: hard acid, oxophilic"),
    ("W", 6):  ElectronicConfig("W",  6, 0, "low", "LMCT",
                                "d⁰: analogous to Mo in enzymes"),
    ("W", 4):  ElectronicConfig("W",  4, 2, "low", "d-d",
                                "d²: strong-field"),
    ("Re", 5): ElectronicConfig("Re", 5, 2, "low", "d-d",
                                "d²: ReO₃-based structures"),
    ("Re", 7): ElectronicConfig("Re", 7, 0, "low", "LMCT",
                                "d⁰: perrhenate, nuclear medicine"),
    ("Os", 2): ElectronicConfig("Os", 2, 6, "low", "d-d / MLCT",
                                "d⁶ LS: strong-field, bipyridyl complexes"),
    ("Os", 4): ElectronicConfig("Os", 4, 4, "low", "d-d",
                                "d⁴: strong-field"),
    ("Ir", 3): ElectronicConfig("Ir", 3, 6, "low", "d-d",
                                "d⁶ LS: extremely inert, octahedral"),
    ("Ir", 4): ElectronicConfig("Ir", 4, 5, "low", "d-d",
                                "d⁵: one unpaired (LS)"),
    ("Pt", 2): ElectronicConfig("Pt", 2, 8, "low", "d-d / MLCT",
                                "d⁸: strong preference for square planar"),
    ("Pt", 4): ElectronicConfig("Pt", 4, 6, "low", "d-d",
                                "d⁶ LS: octahedral, kinetically inert"),
    ("Au", 1): ElectronicConfig("Au", 1, 10, "low", "none (d¹⁰)",
                                "d¹⁰: linear, aurophilic interactions"),
    ("Au", 3): ElectronicConfig("Au", 3, 8, "low", "d-d / LMCT",
                                "d⁸: square planar, strong oxidant"),
    ("Hg", 2): ElectronicConfig("Hg", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰: linear/tetrahedral, spectroscopically silent"),

    # ══════════════════════════════════════════════════════════════════
    # Post-transition metals
    # ══════════════════════════════════════════════════════════════════
    ("Al", 3): ElectronicConfig("Al", 3, 0, "low", "none",
                                "noble gas core, hard Lewis acid", block="s"),
    ("Ga", 3): ElectronicConfig("Ga", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰: Zn-like coordination, Fe³⁺ mimic"),
    ("In", 3): ElectronicConfig("In", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰: soft Lewis acid, flexible CN"),
    ("Tl", 1): ElectronicConfig("Tl", 1, 10, "low", "none (d¹⁰)",
                                "d¹⁰ + 6s²: inert pair, K⁺ mimic"),
    ("Tl", 3): ElectronicConfig("Tl", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰: strong oxidant"),
    ("Sn", 2): ElectronicConfig("Sn", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰ + 5s²: inert pair, hemidirected"),
    ("Sn", 4): ElectronicConfig("Sn", 4, 10, "low", "none (d¹⁰)",
                                "d¹⁰: tetrahedral/octahedral"),
    ("Pb", 2): ElectronicConfig("Pb", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰ + 6s²: inert pair, hemidirected"),
    ("Pb", 4): ElectronicConfig("Pb", 4, 10, "low", "none (d¹⁰)",
                                "d¹⁰: strong oxidant"),
    ("Bi", 3): ElectronicConfig("Bi", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰ + 6s²: inert pair, large CN"),
    ("Sb", 3): ElectronicConfig("Sb", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰ + 5s²: inert pair effect"),
    ("Sb", 5): ElectronicConfig("Sb", 5, 10, "low", "none (d¹⁰)",
                                "d¹⁰: octahedral preference"),
    # Metalloids (appear as HETATM in PDB)
    ("Si", 4): ElectronicConfig("Si", 4, 0, "low", "none",
                                "noble gas core, tetrahedral", block="s"),
    ("Ge", 4): ElectronicConfig("Ge", 4, 10, "low", "none (d¹⁰)",
                                "d¹⁰: tetrahedral, Si analog"),
    ("As", 3): ElectronicConfig("As", 3, 10, "low", "none (d¹⁰)",
                                "d¹⁰: pyramidal / tetrahedral"),
    ("As", 5): ElectronicConfig("As", 5, 10, "low", "none (d¹⁰)",
                                "d¹⁰: arsenate, phosphate mimic"),
    ("Se", 4): ElectronicConfig("Se", 4, 10, "low", "none (d¹⁰)",
                                "d¹⁰: selenocysteine in enzymes"),
    ("Se", 2): ElectronicConfig("Se", 2, 10, "low", "none (d¹⁰)",
                                "d¹⁰: selenolate donor in proteins"),
    ("Te", 4): ElectronicConfig("Te", 4, 10, "low", "none (d¹⁰)",
                                "d¹⁰: tellurite, larger Se analog"),

    # ══════════════════════════════════════════════════════════════════
    # f-block — lanthanides
    # ══════════════════════════════════════════════════════════════════
    ("La", 3): ElectronicConfig("La", 3, 0, "low", "none",
                                "f⁰ d⁰: diamagnetic, colorless, hard Lewis acid",
                                f_count=0, block="f"),
    ("Ce", 2): ElectronicConfig("Ce", 2, 0, "low", "f-f",
                                "f²: rare divalent, strong reductant",
                                f_count=2, block="f"),
    ("Ce", 3): ElectronicConfig("Ce", 3, 0, "low", "f-f (weak)",
                                "f¹: one unpaired f-electron, paramagnetic",
                                f_count=1, block="f"),
    ("Ce", 4): ElectronicConfig("Ce", 4, 0, "low", "LMCT",
                                "f⁰: strong oxidizer, yellow from LMCT",
                                f_count=0, block="f"),
    ("Pr", 2): ElectronicConfig("Pr", 2, 0, "low", "f-f",
                                "f³: rare divalent",
                                f_count=3, block="f"),
    ("Pr", 3): ElectronicConfig("Pr", 3, 0, "low", "f-f (weak)",
                                "f²: green, two unpaired f-electrons",
                                f_count=2, block="f"),
    ("Nd", 2): ElectronicConfig("Nd", 2, 0, "low", "f-f",
                                "f⁴: rare divalent, strong reductant",
                                f_count=4, block="f"),
    ("Nd", 3): ElectronicConfig("Nd", 3, 0, "low", "f-f",
                                "f³: lilac/pink, sharp f-f absorption bands",
                                f_count=3, block="f"),
    ("Pm", 3): ElectronicConfig("Pm", 3, 0, "low", "f-f",
                                "f⁴: radioactive, no stable isotopes",
                                f_count=4, block="f"),
    ("Sm", 3): ElectronicConfig("Sm", 3, 0, "low", "f-f",
                                "f⁵: pale yellow, moderate paramagnet",
                                f_count=5, block="f"),
    ("Sm", 2): ElectronicConfig("Sm", 2, 0, "low", "f-f / 5d←4f",
                                "f⁶: blood-red, strong reductant",
                                f_count=6, block="f"),
    ("Eu", 3): ElectronicConfig("Eu", 3, 0, "low", "f-f",
                                "f⁶: pale pink, strong red luminescence",
                                f_count=6, block="f"),
    ("Eu", 2): ElectronicConfig("Eu", 2, 0, "low", "f-f / 5d←4f",
                                "f⁷: half-filled f-shell, stable divalent Ln",
                                f_count=7, block="f"),
    ("Gd", 2): ElectronicConfig("Gd", 2, 0, "high", "f-f",
                                "f⁸: rare divalent",
                                f_count=8, block="f"),
    ("Gd", 3): ElectronicConfig("Gd", 3, 0, "high", "none (f-f forbidden)",
                                "f⁷: half-filled, highest paramagnetism, MRI contrast",
                                f_count=7, block="f"),
    ("Tb", 3): ElectronicConfig("Tb", 3, 0, "low", "f-f",
                                "f⁸: green luminescence",
                                f_count=8, block="f"),
    ("Tb", 4): ElectronicConfig("Tb", 4, 0, "low", "f-f / LMCT",
                                "f⁷: half-filled, strong oxidant",
                                f_count=7, block="f"),
    ("Dy", 2): ElectronicConfig("Dy", 2, 0, "low", "f-f",
                                "f¹⁰: rare divalent",
                                f_count=10, block="f"),
    ("Dy", 3): ElectronicConfig("Dy", 3, 0, "low", "f-f",
                                "f⁹: pale yellow, highest magnetic moment",
                                f_count=9, block="f"),
    ("Ho", 3): ElectronicConfig("Ho", 3, 0, "low", "f-f",
                                "f¹⁰: yellow, strong paramagnet",
                                f_count=10, block="f"),
    ("Er", 3): ElectronicConfig("Er", 3, 0, "low", "f-f",
                                "f¹¹: pink, NIR luminescence at 1.5 μm",
                                f_count=11, block="f"),
    ("Tm", 2): ElectronicConfig("Tm", 2, 0, "low", "f-f",
                                "f¹³: rare divalent",
                                f_count=13, block="f"),
    ("Tm", 3): ElectronicConfig("Tm", 3, 0, "low", "f-f",
                                "f¹²: pale green",
                                f_count=12, block="f"),
    ("Yb", 3): ElectronicConfig("Yb", 3, 0, "low", "f-f",
                                "f¹³: colorless, NIR luminescence",
                                f_count=13, block="f"),
    ("Yb", 2): ElectronicConfig("Yb", 2, 0, "low", "f-f / 5d←4f",
                                "f¹⁴: diamagnetic divalent, green",
                                f_count=14, block="f"),
    ("Lu", 3): ElectronicConfig("Lu", 3, 0, "low", "none",
                                "f¹⁴ d⁰: diamagnetic, smallest Ln³⁺ radius",
                                f_count=14, block="f"),

    # ══════════════════════════════════════════════════════════════════
    # f-block — actinides
    # ══════════════════════════════════════════════════════════════════
    ("Th", 4): ElectronicConfig("Th", 4, 0, "low", "none",
                                "f⁰ d⁰: diamagnetic, large CN (8-12)",
                                f_count=0, block="f"),
    ("Pa", 4): ElectronicConfig("Pa", 4, 0, "low", "f-f",
                                "f¹: paramagnetic",
                                f_count=1, block="f"),
    ("Pa", 5): ElectronicConfig("Pa", 5, 0, "low", "LMCT",
                                "f⁰: colorless protactinate",
                                f_count=0, block="f"),
    ("U", 4):  ElectronicConfig("U",  4, 0, "low", "f-f",
                                "f²: green, paramagnetic",
                                f_count=2, block="f"),
    ("U", 6):  ElectronicConfig("U",  6, 0, "low", "LMCT",
                                "f⁰: uranyl (UO₂²⁺), yellow fluorescence",
                                f_count=0, block="f"),
    ("Np", 4): ElectronicConfig("Np", 4, 0, "low", "f-f",
                                "f³: paramagnetic",
                                f_count=3, block="f"),
    ("Np", 5): ElectronicConfig("Np", 5, 0, "low", "f-f",
                                "f²: neptunyl (NpO₂⁺)",
                                f_count=2, block="f"),
    ("Pu", 3): ElectronicConfig("Pu", 3, 0, "low", "f-f",
                                "f⁵: blue-violet",
                                f_count=5, block="f"),
    ("Pu", 4): ElectronicConfig("Pu", 4, 0, "low", "f-f",
                                "f⁴: tan/brown",
                                f_count=4, block="f"),
    ("Am", 3): ElectronicConfig("Am", 3, 0, "low", "f-f",
                                "f⁶: pink, luminescent",
                                f_count=6, block="f"),
}


# ─────────────────────────────────────────────────────────────────────
# d-orbital splitting patterns by geometry
# Each entry: (orbital_name, relative_energy, degeneracy)
# ─────────────────────────────────────────────────────────────────────

SPLITTING_PATTERNS: dict[str, list[tuple[str, float, int]]] = {
    "octahedral": [
        ("t₂g (dxy, dxz, dyz)", 0.0, 3),
        ("eg (dz², dx²-y²)", 1.0, 2),
    ],
    "tetrahedral": [
        ("e (dz², dx²-y²)", 0.0, 2),
        ("t₂ (dxy, dxz, dyz)", 0.44, 3),  # Δt ≈ 4/9 Δo
    ],
    "square planar": [
        ("dxz, dyz", 0.0, 2),
        ("dz²", 0.3, 1),
        ("dxy", 0.7, 1),
        ("dx²-y²", 1.2, 1),
    ],
    "square pyramidal": [
        ("dxz, dyz", 0.0, 2),
        ("dxy", 0.35, 1),
        ("dz²", 0.65, 1),
        ("dx²-y²", 1.0, 1),
    ],
    "trigonal bipyramidal": [
        ("dz²", 0.0, 1),
        ("dx²-y², dxy", 0.3, 2),
        ("dxz, dyz", 0.7, 2),
    ],
    # ── Distorted geometries — degeneracies partially lifted ──
    # These MUST differ from ideal to make comparisons meaningful.
    "distorted tetrahedral": [
        # e → split into dz² / dx²-y², t₂ → split into dxy / (dxz,dyz)
        ("dz²", 0.0, 1),
        ("dx²-y²", 0.08, 1),
        ("dxy", 0.40, 1),
        ("dxz, dyz", 0.50, 2),
    ],
    "distorted octahedral": [
        # Tetragonal distortion (Jahn-Teller or mixed donors):
        # t₂g splits into dxy (lower) + dxz,dyz (upper)
        # eg splits into dz² (lower) + dx²-y² (upper)
        ("dxy", 0.0, 1),
        ("dxz, dyz", 0.08, 2),
        ("dz²", 0.92, 1),
        ("dx²-y²", 1.05, 1),
    ],
    "irregular octahedral": [
        # More severe distortion — all levels split
        ("dxy", 0.0, 1),
        ("dxz", 0.06, 1),
        ("dyz", 0.12, 1),
        ("dz²", 0.88, 1),
        ("dx²-y²", 1.05, 1),
    ],
    "distorted square planar": [
        ("dxz", 0.0, 1),
        ("dyz", 0.05, 1),
        ("dz²", 0.28, 1),
        ("dxy", 0.68, 1),
        ("dx²-y²", 1.15, 1),
    ],
    "distorted trigonal bipyramidal": [
        ("dz²", 0.0, 1),
        ("dx²-y²", 0.25, 1),
        ("dxy", 0.35, 1),
        ("dxz", 0.65, 1),
        ("dyz", 0.75, 1),
    ],
    "distorted square pyramidal": [
        ("dxz", 0.0, 1),
        ("dyz", 0.06, 1),
        ("dxy", 0.32, 1),
        ("dz²", 0.60, 1),
        ("dx²-y²", 1.0, 1),
    ],
    "seesaw / intermediate": [
        ("dxz", 0.0, 1),
        ("dyz", 0.10, 1),
        ("dz²", 0.35, 1),
        ("dxy", 0.65, 1),
        ("dx²-y²", 1.0, 1),
    ],
    "trigonal prismatic": [
        ("dz²", 0.0, 1),
        ("dxz, dyz", 0.25, 2),
        ("dxy, dx²-y²", 0.8, 2),
    ],
    "linear": [
        ("dxy, dx²-y²", 0.0, 2),
        ("dxz, dyz", 0.2, 2),
        ("dz²", 1.0, 1),
    ],
}


# ─────────────────────────────────────────────────────────────────────
# f-orbital splitting patterns by geometry
# Each entry: (orbital_name, relative_energy, degeneracy)
# NOTE: f-orbital splitting in Ln³⁺ is ~100× smaller than d-orbital
# splitting.  The energy scale is relative within each diagram.
# ─────────────────────────────────────────────────────────────────────

F_SPLITTING_PATTERNS: dict[str, list[tuple[str, float, int]]] = {
    "octahedral": [
        ("a₂u (fxyz)", 0.0, 1),
        ("t₂u (fx³, fy³, fz³)", 0.25, 3),
        ("t₁u (fx(y²−z²), fy(z²−x²), fz(x²−y²))", 1.0, 3),
    ],
    "square antiprismatic": [
        # Common CN=8 for lanthanides
        ("a₂", 0.0, 1),
        ("e₁", 0.15, 2),
        ("b₁ + b₂", 0.4, 2),
        ("e₂", 0.8, 2),
    ],
    "tricapped trigonal prismatic": [
        # Common CN=9 for light lanthanides
        ("a₁'", 0.0, 1),
        ("e'", 0.2, 2),
        ("e''", 0.5, 2),
        ("a₁'' + a₂''", 0.9, 2),
    ],
    "tetrahedral": [
        ("a₂ (fxyz)", 0.0, 1),
        ("t₁ (fx³, fy³, fz³)", 0.4, 3),
        ("t₂ (fx(y²−z²), fy(z²−x²), fz(x²−y²))", 0.7, 3),
    ],
    "cubic": [
        # CN=8 cubic
        ("a₂u (fxyz)", 0.0, 1),
        ("t₂u (fx³, fy³, fz³)", 0.5, 3),
        ("t₁u", 1.0, 3),
    ],
    # Fallback for unrecognized geometries
    "default": [
        ("lower f-levels", 0.0, 4),
        ("upper f-levels", 0.6, 3),
    ],
    # CN-based aliases (common for lanthanides)
    "8-coordinate": [
        # Approximate square antiprismatic
        ("a₂", 0.0, 1),
        ("e₁", 0.15, 2),
        ("b₁ + b₂", 0.4, 2),
        ("e₂", 0.8, 2),
    ],
    "9-coordinate": [
        # Approximate tricapped trigonal prismatic
        ("a₁'", 0.0, 1),
        ("e'", 0.2, 2),
        ("e''", 0.5, 2),
        ("a₁'' + a₂''", 0.9, 2),
    ],
    "7-coordinate": [
        # Approximate capped octahedral / pentagonal bipyramidal
        ("a₁", 0.0, 1),
        ("e₁", 0.15, 2),
        ("e₂", 0.45, 2),
        ("e₃", 0.85, 2),
    ],
    "10-coordinate": [
        # Approximate bicapped square antiprismatic
        ("a₁", 0.0, 1),
        ("e₁", 0.1, 2),
        ("b₁", 0.3, 1),
        ("e₂", 0.55, 2),
        ("a₂", 0.85, 1),
    ],
    "6-coordinate": [
        # Same as octahedral
        ("a₂u (fxyz)", 0.0, 1),
        ("t₂u (fx³, fy³, fz³)", 0.25, 3),
        ("t₁u (fx(y²−z²), fy(z²−x²), fz(x²−y²))", 1.0, 3),
    ],
}


# ─────────────────────────────────────────────────────────────────────
# API functions
# ─────────────────────────────────────────────────────────────────────

def get_config(element: str, oxidation_state: int) -> Optional[ElectronicConfig]:
    """Look up electronic configuration for a metal ion.

    If the exact (element, oxidation_state) pair is not in CONFIGS,
    tries the most common oxidation state for f-block elements (+3).
    """
    config = CONFIGS.get((element, oxidation_state))
    if config:
        return config

    # Fallback: for f-block (lanthanides/actinides), +3 is almost universal
    LANTHANIDES = {"La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu"}
    ACTINIDES = {"Th","Pa","U","Np","Pu","Am"}
    if element in LANTHANIDES or element in ACTINIDES:
        # Try +3 first, then scan for any available config
        for try_ox in [3, 4, 2, 5, 6]:
            fallback = CONFIGS.get((element, try_ox))
            if fallback:
                return fallback

    return None


def get_splitting(geometry: str) -> list[tuple[str, float, int]]:
    """
    Get d-orbital splitting pattern for a geometry.
    Falls back to the closest match if exact geometry not found.
    """
    geo_lower = geometry.lower()

    if geo_lower in SPLITTING_PATTERNS:
        return SPLITTING_PATTERNS[geo_lower]

    for key, pattern in SPLITTING_PATTERNS.items():
        if key in geo_lower or geo_lower in key:
            return pattern

    return SPLITTING_PATTERNS["octahedral"]


def get_f_splitting(geometry: str) -> list[tuple[str, float, int]]:
    """
    Get f-orbital splitting pattern for a geometry.
    Falls back to a generic two-level split if not found.
    """
    geo_lower = geometry.lower()

    if geo_lower in F_SPLITTING_PATTERNS:
        return F_SPLITTING_PATTERNS[geo_lower]

    for key, pattern in F_SPLITTING_PATTERNS.items():
        if key in geo_lower or geo_lower in key:
            return pattern

    # Default: most lanthanide CN≥7, use a generic split
    return F_SPLITTING_PATTERNS["default"]


def fill_orbitals(
    electron_count: int,
    splitting: list[tuple[str, float, int]],
    high_spin: bool = False,
) -> list[tuple[str, float, int, int]]:
    """
    Fill electrons into orbital levels.

    Returns list of (label, energy, degeneracy, n_electrons).
    Applies aufbau + Hund for high spin, or fill-lowest for low spin.
    """
    if high_spin:
        filled = []
        remaining = electron_count
        # Pass 1: one per orbital
        for label, energy, deg in splitting:
            n = min(deg, remaining)
            filled.append((label, energy, deg, n))
            remaining -= n
        # Pass 2: pair up
        if remaining > 0:
            new_filled = []
            for label, energy, deg, n_e in filled:
                add = min(deg - n_e, remaining) if remaining > 0 else 0
                new_filled.append((label, energy, deg, n_e + add))
                remaining -= add
            filled = new_filled
        return filled
    else:
        filled = []
        remaining = electron_count
        for label, energy, deg in splitting:
            n = min(2 * deg, remaining)
            filled.append((label, energy, deg, n))
            remaining -= n
        return filled

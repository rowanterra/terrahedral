"""
Catalytic analysis for coordination sites.

Provides:
  - Catalytic function inference from donor set + geometry + metal
  - Entatic state quantification (distortion from thermodynamic preference)
  - CatalyticCycle model for multi-step enzyme mechanisms

⚠ DISCLAIMER: Catalytic function predictions are heuristic inferences
based on published coordination motifs and are NOT derived from the CIF
structural data.  They should be treated as hypotheses for further
investigation, not definitive assignments.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from terrahedral.core import CoordinationSite


# ─────────────────────────────────────────────────────────────────────
# Catalytic function inference
# ─────────────────────────────────────────────────────────────────────

# Known coordination motifs → probable catalytic function.
# Keys: (element, frozenset of donor elements, CN range)
# Each entry: list of (function, confidence, description)

@dataclass
class FunctionPrediction:
    """A predicted catalytic function with confidence level."""
    function: str           # e.g. "oxidase", "hydrolase"
    confidence: str         # "high", "moderate", "low"
    basis: str              # why we think this
    examples: str = ""      # known enzymes with this motif
    disclaimer: str = (
        "Inferred from coordination motif — not derived from CIF data. "
        "Treat as hypothesis for further investigation."
    )


# Motif database: maps (element, donor_elements) → predictions
_MOTIF_DB: list[dict] = [
    # ── Iron ──
    {
        "element": "Fe",
        "donors": {"N", "N", "O"},  # will match subsets
        "donor_formula": "N2-3O1-3",
        "cn": (5, 6),
        "predictions": [
            FunctionPrediction(
                "2-oxoglutarate-dependent oxygenase",
                "high",
                "2-His-1-carboxylate facial triad (HX...H...E/D motif) with Fe²⁺",
                "AsqJ, TauD, FTO, AlkB, PHD2, P4H",
            ),
            FunctionPrediction(
                "halogenase",
                "moderate",
                "Same facial triad but with halide (Cl⁻) replacing one water",
                "SyrB2, WelO5, BesD",
            ),
        ],
    },
    {
        "element": "Fe",
        "donors": {"N", "S"},
        "donor_formula": "N3-4S1",
        "cn": (4, 5, 6),
        "predictions": [
            FunctionPrediction(
                "nitrile hydratase",
                "moderate",
                "Fe³⁺ with Cys + His donors, non-heme non-corrin",
                "NHase",
            ),
        ],
    },
    {
        "element": "Fe",
        "donors": {"S"},
        "donor_formula": "S4",
        "cn": (4,),
        "predictions": [
            FunctionPrediction(
                "electron transfer (rubredoxin-type)",
                "high",
                "Tetrahedral Fe-S₄ from 4 Cys residues",
                "Rubredoxin, desulforedoxin",
            ),
        ],
    },
    # ── Copper ──
    {
        "element": "Cu",
        "donors": {"N", "S"},
        "donor_formula": "N2S1 or N2S1O1",
        "cn": (3, 4),
        "predictions": [
            FunctionPrediction(
                "electron transfer (type 1 / blue copper)",
                "high",
                "Trigonal/tetrahedral Cu with 2 His + Cys(Sγ) ± Met(Sδ)",
                "Plastocyanin, azurin, stellacyanin, laccase T1",
            ),
        ],
    },
    {
        "element": "Cu",
        "donors": {"N"},
        "donor_formula": "N3-4",
        "cn": (4, 5),
        "predictions": [
            FunctionPrediction(
                "oxidase / O₂ activation (type 2 copper)",
                "moderate",
                "Square planar / pyramidal Cu²⁺ with His donors",
                "Galactose oxidase, amine oxidase, Cu/Zn SOD",
            ),
        ],
    },
    # ── Zinc ──
    {
        "element": "Zn",
        "donors": {"N", "O"},
        "donor_formula": "N2-3O1-2",
        "cn": (4, 5),
        "predictions": [
            FunctionPrediction(
                "hydrolase (Lewis acid catalysis)",
                "high",
                "Tetrahedral/5-coord Zn²⁺ with His + Asp/Glu + water",
                "Carbonic anhydrase, carboxypeptidase, thermolysin, MMP",
            ),
        ],
    },
    {
        "element": "Zn",
        "donors": {"S"},
        "donor_formula": "S3-4",
        "cn": (4,),
        "predictions": [
            FunctionPrediction(
                "structural / DNA binding (zinc finger)",
                "high",
                "Tetrahedral Zn²⁺ with Cys₄ or Cys₂His₂",
                "Zinc finger proteins, GATA, nuclear receptors",
            ),
        ],
    },
    # ── Manganese ──
    {
        "element": "Mn",
        "donors": {"N", "O"},
        "donor_formula": "N1-2O3-5",
        "cn": (5, 6),
        "predictions": [
            FunctionPrediction(
                "superoxide dismutase",
                "moderate",
                "Mn²⁺/³⁺ with His + Asp in trigonal bipyramidal/octahedral",
                "MnSOD",
            ),
            FunctionPrediction(
                "water oxidation (OEC-like)",
                "low",
                "High-valent Mn cluster with bridging oxo",
                "Photosystem II OEC (Mn₄CaO₅)",
            ),
        ],
    },
    # ── Nickel ──
    {
        "element": "Ni",
        "donors": {"S", "N"},
        "donor_formula": "N2S2 or S4",
        "cn": (4, 5, 6),
        "predictions": [
            FunctionPrediction(
                "hydrogenase / CO dehydrogenase",
                "moderate",
                "Ni with Cys thiolates, often square planar",
                "[NiFe]-hydrogenase, CODH/ACS",
            ),
        ],
    },
    # ── Cobalt ──
    {
        "element": "Co",
        "donors": {"N"},
        "donor_formula": "N4-5",
        "cn": (5, 6),
        "predictions": [
            FunctionPrediction(
                "isomerase / methyl transfer (B₁₂-dependent)",
                "moderate",
                "Octahedral Co³⁺ with equatorial N₄ (corrin-like)",
                "Methylmalonyl-CoA mutase, methionine synthase",
            ),
        ],
    },
    # ── Molybdenum ──
    {
        "element": "Mo",
        "donors": {"S", "O"},
        "donor_formula": "S2O1-2 or S1O2-3",
        "cn": (4, 5, 6),
        "predictions": [
            FunctionPrediction(
                "oxo-transfer (molybdopterin enzyme)",
                "high",
                "Mo with dithiolene (S₂) + oxo/hydroxo + protein ligand",
                "Sulfite oxidase, xanthine oxidase, DMSO reductase, nitrate reductase",
            ),
        ],
    },
    # ── Vanadium ──
    {
        "element": "V",
        "donors": {"N", "O"},
        "donor_formula": "N1-3O2-4",
        "cn": (5, 6),
        "predictions": [
            FunctionPrediction(
                "haloperoxidase",
                "moderate",
                "Trigonal bipyramidal V⁵⁺ with His + oxo groups",
                "V-chloroperoxidase, V-bromoperoxidase",
            ),
        ],
    },
    # ── Tungsten ──
    {
        "element": "W",
        "donors": {"S", "O"},
        "donor_formula": "S2-4O1-2",
        "cn": (5, 6),
        "predictions": [
            FunctionPrediction(
                "oxo-transfer (tungstoenzyme)",
                "moderate",
                "W with dithiolene(s) + oxo, analogous to Mo enzymes",
                "Aldehyde oxidoreductase (AOR), formate dehydrogenase",
            ),
        ],
    },
    # ── Lanthanides ──
    {
        "element": "Nd",
        "donors": {"O", "N"},
        "donor_formula": "O5-9N0-3",
        "cn": (7, 8, 9),
        "predictions": [
            FunctionPrediction(
                "methanol dehydrogenase (XoxF-type)",
                "moderate",
                "Ln³⁺ in PQQ-dependent alcohol oxidation, replacing Ca²⁺",
                "XoxF-MDH, ExaF",
            ),
            FunctionPrediction(
                "lanthanide-binding / storage (lanmodulin)",
                "high",
                "EF-hand-like Ln³⁺ binding with carboxylate-rich coordination",
                "Lanmodulin (LanM), Hans-LanM",
            ),
        ],
    },
]

# Copy Nd entry for all lanthanides
for _ln in ["La", "Ce", "Pr", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu"]:
    _motif = dict(_MOTIF_DB[-1])
    _motif["element"] = _ln
    _MOTIF_DB.append(_motif)


def infer_function(site: CoordinationSite) -> list[FunctionPrediction]:
    """
    Predict catalytic function from the coordination site.

    Returns a list of FunctionPrediction objects, sorted by confidence.

    ⚠ These are heuristic inferences, NOT experimental assignments.
    """
    el = site.metal.element
    cn = site.coordination_number
    donor_els = {lig.element for lig in site.ligands}

    results = []
    for motif in _MOTIF_DB:
        if motif["element"] != el:
            continue
        if cn not in motif["cn"]:
            continue
        # Check donor overlap — motif donors should be a subset of actual donors
        if not motif["donors"].issubset(donor_els):
            continue
        results.extend(motif["predictions"])

    # Sort by confidence
    order = {"high": 0, "moderate": 1, "low": 2}
    results.sort(key=lambda p: order.get(p.confidence, 3))
    return results


# ─────────────────────────────────────────────────────────────────────
# Entatic state analysis
# ─────────────────────────────────────────────────────────────────────

# Preferred (thermodynamic) geometry for metal ions in the absence of
# protein constraints.  These are the "free ion" preferences.
PREFERRED_GEOMETRY: dict[tuple[str, int], str] = {
    # Cu
    ("Cu", 1): "linear",            # d¹⁰ → linear or trigonal
    ("Cu", 2): "square planar",     # d⁹ → Jahn-Teller elongated oct or sq pl
    # Fe
    ("Fe", 2): "octahedral",        # d⁶ HS → regular octahedral
    ("Fe", 3): "octahedral",        # d⁵ HS → regular octahedral
    # Zn
    ("Zn", 2): "tetrahedral",       # d¹⁰ → tetrahedral
    # Mn
    ("Mn", 2): "octahedral",
    ("Mn", 3): "octahedral",
    # Co
    ("Co", 2): "octahedral",        # but tetrahedral also common
    ("Co", 3): "octahedral",
    # Ni
    ("Ni", 2): "square planar",     # in strong field; octahedral in weak
    # Pt
    ("Pt", 2): "square planar",     # strong d⁸ preference
    ("Pt", 4): "octahedral",
    # Mo
    ("Mo", 6): "octahedral",
    ("Mo", 4): "octahedral",
    # V
    ("V", 5): "trigonal bipyramidal",
    ("V", 4): "square pyramidal",
}


@dataclass
class EntaicAnalysis:
    """Entatic state / rack mechanism analysis."""
    actual_geometry: str
    preferred_geometry: str
    is_entatic: bool        # True if actual ≠ preferred
    tau_actual: Optional[float]
    tau_ideal: Optional[float]
    distortion_description: str
    catalytic_significance: str
    disclaimer: str = (
        "Entatic state inference is based on comparison of crystallographic "
        "geometry to thermodynamic preference for the free metal ion. "
        "Actual protein-induced strain depends on many factors not captured here."
    )


def analyze_entatic(
    site: CoordinationSite,
    oxidation_state: int = 2,
) -> Optional[EntaicAnalysis]:
    """
    Analyze whether the coordination site shows entatic state features.

    The entatic state / rack mechanism (Vallee & Williams 1968, Malmström 1965)
    proposes that proteins distort metal coordination away from thermodynamic
    preference to optimize catalytic function.

    Classic example: blue copper proteins force Cu²⁺ into trigonal/
    tetrahedral geometry (preferred by Cu⁺) rather than square planar
    (preferred by Cu²⁺), facilitating fast electron transfer by
    minimizing reorganization energy.
    """
    from terrahedral.analysis import tau4, tau5, classify, ideal_geometry

    el = site.metal.element
    preferred = PREFERRED_GEOMETRY.get((el, oxidation_state))
    if preferred is None:
        return None

    actual = site.geometry or classify(site)
    is_entatic = actual.lower().replace("distorted ", "") != preferred.lower()

    # Compute tau values
    cn = site.coordination_number
    t_actual = tau4(site) if cn == 4 else (tau5(site) if cn == 5 else None)

    # Describe the distortion
    if not is_entatic:
        desc = f"Geometry matches thermodynamic preference ({preferred}). No significant entatic strain."
        signif = "Standard coordination — reorganization energy not minimized by protein."
    else:
        desc = (
            f"Protein enforces {actual} geometry on {el}{oxidation_state}+, "
            f"which thermodynamically prefers {preferred}."
        )

        # Specific catalytic significance
        if el == "Cu" and oxidation_state == 2 and "tetra" in actual.lower():
            signif = (
                "Classic blue copper entatic state: Cu²⁺ held in Cu⁺-like "
                "tetrahedral geometry minimizes reorganization energy for "
                "fast electron transfer (Marcus theory)."
            )
        elif el == "Fe" and "trigonal" in actual.lower():
            signif = (
                "Distortion toward trigonal bipyramidal may create an open "
                "coordination site for substrate/O₂ binding."
            )
        elif el == "Fe" and "square" in actual.lower() and cn == 5:
            signif = (
                "Square pyramidal geometry with open axial site is typical "
                "of the 2-His-1-carboxylate triad poised for O₂ activation."
            )
        elif cn != {"tetrahedral": 4, "square planar": 4, "octahedral": 6,
                     "trigonal bipyramidal": 5, "square pyramidal": 5,
                     "linear": 2}.get(preferred, cn):
            signif = (
                f"Change in coordination number ({cn} vs typical "
                f"{preferred} CN) suggests protein-enforced geometry "
                f"to tune redox potential or substrate access."
            )
        else:
            signif = (
                f"Distortion from {preferred} suggests protein-tuned "
                f"geometry for optimized catalytic function."
            )

    return EntaicAnalysis(
        actual_geometry=actual,
        preferred_geometry=preferred,
        is_entatic=is_entatic,
        tau_actual=t_actual,
        tau_ideal=None,
        distortion_description=desc,
        catalytic_significance=signif,
    )


# ─────────────────────────────────────────────────────────────────────
# Catalytic cycle model
# ─────────────────────────────────────────────────────────────────────

@dataclass
class CycleStep:
    """One step in a catalytic cycle."""
    name: str                       # e.g. "Resting state", "Ferryl Fe(IV)"
    site: CoordinationSite          # coordination at this step
    oxidation_state: int = 2
    spin: int = 0                   # total spin quantum number S
    description: str = ""           # e.g. "2-His-1-carboxylate + 3 waters"
    pdb_id: str = ""                # source structure if known
    is_inferred: bool = False       # True if not from a crystal structure


@dataclass
class CatalyticCycle:
    """
    A multi-step catalytic mechanism with coordination geometry at each step.

    Build from CIF files (one per step) or programmatically.
    """
    title: str
    steps: list[CycleStep] = field(default_factory=list)
    reference: str = ""             # e.g. "Guo / Chang · JACS 2020"
    pdb_ids: str = ""               # e.g. "5DAO → 5DB0"

    def add_step(
        self,
        name: str,
        site: CoordinationSite,
        oxidation_state: int = 2,
        spin: int = 0,
        description: str = "",
        pdb_id: str = "",
        is_inferred: bool = False,
    ) -> None:
        """Add a step to the cycle."""
        site.classify_geometry()
        self.steps.append(CycleStep(
            name=name,
            site=site,
            oxidation_state=oxidation_state,
            spin=spin,
            description=description,
            pdb_id=pdb_id,
            is_inferred=is_inferred,
        ))

    def render_html(self) -> str:
        """Render the catalytic cycle as an interactive HTML widget."""
        from terrahedral.renderers.html import render_cycle
        return render_cycle(self)

    @property
    def metal_element(self) -> str:
        if self.steps:
            return self.steps[0].site.metal.element
        return "?"

"""
Catalytic function prediction from metal coordination geometry.

Maps (metal, oxidation state, coordination geometry, donor set) → probable
catalytic function, known enzyme families, and mineral analogs.

The core idea: metalloenzymes and minerals that share coordination motifs
and Eh-pH stability windows often perform analogous redox chemistry.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

from terrahedral.core import CoordinationSite


@dataclass
class CatalyticPrediction:
    """Predicted catalytic function for a coordination site."""
    function: str                    # e.g. "electron transfer", "O₂ activation"
    confidence: str = "moderate"     # "high", "moderate", "speculative"
    mechanism: str = ""              # e.g. "inner-sphere ET", "radical rebound"
    enzyme_families: list[str] = field(default_factory=list)
    mineral_analogs: list[str] = field(default_factory=list)
    eh_window: tuple[float, float] = (0.0, 0.0)  # (Eh_min, Eh_max) in V vs SHE
    ph_window: tuple[float, float] = (0.0, 14.0)  # (pH_min, pH_max)
    notes: str = ""


# ─────────────────────────────────────────────────────────────────────
# Knowledge base: (element, oxidation_states, geometry_pattern, donor_pattern) → function
# ─────────────────────────────────────────────────────────────────────

_RULES: list[dict] = [
    # ── Iron ──
    {
        "elements": ["Fe"],
        "ox_states": [2, 3],
        "geometries": ["octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: "N" in ds and "O" in ds,
        "prediction": CatalyticPrediction(
            function="O₂ activation / hydroxylation",
            confidence="high",
            mechanism="Radical rebound: FeII + O₂ → FeIV=O → H-atom abstraction → OH rebound",
            enzyme_families=["Non-heme Fe dioxygenases", "α-ketoglutarate-dependent hydroxylases",
                             "Rieske dioxygenases", "TauD"],
            mineral_analogs=["Goethite (α-FeOOH)", "Lepidocrocite (γ-FeOOH)"],
            eh_window=(-0.1, 0.8),
            ph_window=(4.0, 9.0),
            notes="FeII/FeIII cycling with O₂. N/O donor set typical of 2-His-1-carboxylate facial triad.",
        ),
    },
    {
        "elements": ["Fe"],
        "ox_states": [2, 3],
        "geometries": ["tetrahedral", "distorted tetrahedral"],
        "donor_pattern": lambda ds: ds.count("S") >= 3,
        "prediction": CatalyticPrediction(
            function="Electron transfer",
            confidence="high",
            mechanism="Outer-sphere single-electron transfer via Fe-S cluster",
            enzyme_families=["Ferredoxins", "Rubredoxins", "Fe-S cluster proteins",
                             "Rieske proteins", "Nitrogenase P-cluster"],
            mineral_analogs=["Pyrite (FeS₂)", "Mackinawite (FeS)", "Greigite (Fe₃S₄)"],
            eh_window=(-0.5, 0.1),
            ph_window=(5.0, 9.0),
            notes="Tetrahedral FeS₄ is the universal electron shuttle in biology. "
                  "Mineral analogs in hydrothermal vents may have been primordial ET catalysts.",
        ),
    },
    {
        "elements": ["Fe"],
        "ox_states": [2, 3],
        "geometries": ["square pyramidal", "distorted square pyramidal", "octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: ds.count("N") >= 4,
        "prediction": CatalyticPrediction(
            function="O₂ binding / electron transfer / peroxidase",
            confidence="high",
            mechanism="Axial O₂ binding to FeII in porphyrin plane; reversible or with H₂O₂ activation",
            enzyme_families=["Hemoglobin/Myoglobin", "Cytochrome P450", "Peroxidases",
                             "Cytochrome c", "Catalase"],
            mineral_analogs=["Hematite (Fe₂O₃)", "Magnetite (Fe₃O₄)"],
            eh_window=(-0.2, 1.0),
            ph_window=(5.0, 9.0),
            notes="Porphyrin-coordinated Fe (heme). N₄ equatorial + axial ligand(s). "
                  "The most versatile catalytic platform in biology.",
        ),
    },
    # ── Copper ──
    {
        "elements": ["Cu"],
        "ox_states": [1, 2],
        "geometries": ["distorted tetrahedral", "trigonal pyramidal", "tetrahedral"],
        "donor_pattern": lambda ds: "N" in ds and ("S" in ds or "O" in ds),
        "prediction": CatalyticPrediction(
            function="Electron transfer (type 1 / blue copper)",
            confidence="high",
            mechanism="Outer-sphere ET; entatic state enforces tetrahedral CuII (normally prefers sq. planar)",
            enzyme_families=["Plastocyanin", "Azurin", "Stellacyanin", "Rusticyanin",
                             "Amicyanin", "Laccase (T1 site)"],
            mineral_analogs=["Chalcopyrite (CuFeS₂)", "Covellite (CuS)"],
            eh_window=(0.2, 0.8),
            ph_window=(4.0, 9.0),
            notes="Classic entatic state: protein forces Cu into geometry unfavorable for CuII "
                  "but ideal for rapid ET. High reduction potential (200–800 mV).",
        ),
    },
    {
        "elements": ["Cu"],
        "ox_states": [1, 2],
        "geometries": ["square planar", "distorted square planar", "square pyramidal",
                        "distorted square pyramidal"],
        "donor_pattern": lambda ds: "N" in ds,
        "prediction": CatalyticPrediction(
            function="O₂ activation / oxidase / oxygenase",
            confidence="high",
            mechanism="CuI + O₂ → CuII-superoxo or bis-Cu peroxo; substrate oxidation",
            enzyme_families=["Tyrosinase", "Catechol oxidase", "Amine oxidase",
                             "Galactose oxidase", "Cu-Zn SOD (Cu site)"],
            mineral_analogs=["Cuprite (Cu₂O)", "Tenorite (CuO)", "Malachite"],
            eh_window=(0.0, 0.6),
            ph_window=(5.0, 9.0),
            notes="Type 2/3 Cu. Square planar or pyramidal geometry supports CuII stabilization "
                  "and O₂-derived intermediate chemistry.",
        ),
    },
    # ── Zinc ──
    {
        "elements": ["Zn"],
        "ox_states": [2],
        "geometries": ["tetrahedral", "distorted tetrahedral"],
        "donor_pattern": lambda ds: True,
        "prediction": CatalyticPrediction(
            function="Lewis acid catalysis / hydrolysis",
            confidence="high",
            mechanism="ZnII polarizes substrate (H₂O, CO₂, peptide); no redox. "
                      "Open coordination site for substrate binding.",
            enzyme_families=["Carbonic anhydrase", "Carboxypeptidase", "Thermolysin",
                             "Alcohol dehydrogenase", "Matrix metalloproteinases",
                             "Zinc finger proteins (structural)"],
            mineral_analogs=["Sphalerite (ZnS)", "Smithsonite (ZnCO₃)", "Hydrozincite"],
            eh_window=(-0.5, 0.5),
            ph_window=(5.0, 10.0),
            notes="Zn²⁺ is redox-inert (d¹⁰). Catalysis is purely Lewis acid: "
                  "lowers pKa of bound water, activates electrophiles.",
        ),
    },
    # ── Manganese ──
    {
        "elements": ["Mn"],
        "ox_states": [2, 3, 4],
        "geometries": ["octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: "O" in ds,
        "prediction": CatalyticPrediction(
            function="Water oxidation / O₂ evolution / radical scavenging",
            confidence="high",
            mechanism="Multi-electron oxidation: MnII→MnIII→MnIV. "
                      "OEC: Mn₄CaO₅ cluster cycles through S-states.",
            enzyme_families=["Photosystem II (OEC)", "Mn-SOD", "Mn-catalase",
                             "Arginase", "Oxalate decarboxylase"],
            mineral_analogs=["Pyrolusite (MnO₂)", "Birnessite (δ-MnO₂)",
                             "Rhodochrosite (MnCO₃)", "Manganite (MnOOH)"],
            eh_window=(-0.2, 1.2),
            ph_window=(4.0, 10.0),
            notes="Mn can access +2/+3/+4 oxidation states in biology. "
                  "MnO₂ minerals are natural water oxidation catalysts.",
        ),
    },
    # ── Cobalt ──
    {
        "elements": ["Co"],
        "ox_states": [2, 3],
        "geometries": ["octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: "N" in ds,
        "prediction": CatalyticPrediction(
            function="Radical generation / alkyl transfer / isomerization",
            confidence="high",
            mechanism="Homolytic Co-C bond cleavage → 5'-deoxyadenosyl radical; "
                      "or heterolytic methyl transfer",
            enzyme_families=["B₁₂-dependent enzymes", "Methylmalonyl-CoA mutase",
                             "Methionine synthase", "Diol dehydratase"],
            mineral_analogs=["Cobaltite (CoAsS)", "Linnaeite (Co₃S₄)"],
            eh_window=(-0.2, 0.5),
            ph_window=(5.0, 9.0),
            notes="Corrin ring (modified porphyrin) coordinates Co. "
                  "Only biological system that uses metal-carbon bonds.",
        ),
    },
    # ── Nickel ──
    {
        "elements": ["Ni"],
        "ox_states": [2, 3],
        "geometries": ["square planar", "distorted square planar", "octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: True,
        "prediction": CatalyticPrediction(
            function="Thiol/H₂ metabolism / CO₂ fixation",
            confidence="moderate",
            mechanism="NiII/NiIII cycling; H₂ activation at Ni-Fe center; "
                      "CO₂ reduction to CO at Ni-Fe-S cluster",
            enzyme_families=["[NiFe]-hydrogenase", "CODH/ACS", "Urease",
                             "Ni-SOD", "Methyl-CoM reductase (F430)"],
            mineral_analogs=["Pentlandite ((Fe,Ni)₉S₈)", "Millerite (NiS)", "Awaruite (Ni₃Fe)"],
            eh_window=(-0.6, 0.2),
            ph_window=(5.0, 9.0),
            notes="Ni enzymes are concentrated in anaerobic/ancient metabolisms "
                  "(H₂, CO, CH₄). Strong link to hydrothermal vent mineralogy.",
        ),
    },
    # ── Molybdenum ──
    {
        "elements": ["Mo"],
        "ox_states": [4, 5, 6],
        "geometries": ["square pyramidal", "distorted square pyramidal", "octahedral", "distorted octahedral"],
        "donor_pattern": lambda ds: "S" in ds or "O" in ds,
        "prediction": CatalyticPrediction(
            function="Oxo-transfer / oxygen atom transfer",
            confidence="high",
            mechanism="MoVI=O + substrate → MoIV + substrate-O; "
                      "regenerated by H₂O or O₂",
            enzyme_families=["Sulfite oxidase", "Xanthine oxidase", "DMSO reductase",
                             "Nitrate reductase", "Aldehyde oxidase"],
            mineral_analogs=["Molybdenite (MoS₂)", "Wulfenite (PbMoO₄)"],
            eh_window=(-0.3, 0.6),
            ph_window=(5.0, 9.0),
            notes="Mo uses the molybdopterin cofactor (pyranopterin dithiolene). "
                  "2-electron oxo-transfer is the signature reaction.",
        ),
    },
    # ── Vanadium ──
    {
        "elements": ["V"],
        "ox_states": [3, 4, 5],
        "geometries": ["octahedral", "distorted octahedral", "square pyramidal", "trigonal bipyramidal"],
        "donor_pattern": lambda ds: True,
        "prediction": CatalyticPrediction(
            function="Haloperoxidase / nitrogenase / phosphatase mimic",
            confidence="moderate",
            mechanism="VV-peroxo intermediate halogenates substrates; "
                      "V-nitrogenase reduces N₂ via V-Fe-S cluster",
            enzyme_families=["V-haloperoxidase", "V-nitrogenase", "V-chloroperoxidase"],
            mineral_analogs=["Vanadinite (Pb₅(VO₄)₃Cl)", "Carnotite (K₂(UO₂)₂(VO₄)₂)"],
            eh_window=(-0.3, 0.8),
            ph_window=(3.0, 9.0),
            notes="V is a phosphate mimic (vanadate ≈ phosphate geometry). "
                  "Concentrated in tunicates and Amanita mushrooms.",
        ),
    },
    # ── Platinum ──
    {
        "elements": ["Pt"],
        "ox_states": [2, 4],
        "geometries": ["square planar", "distorted square planar"],
        "donor_pattern": lambda ds: True,
        "prediction": CatalyticPrediction(
            function="DNA cross-linking / catalytic hydrogenation (synthetic)",
            confidence="high",
            mechanism="PtII square planar → aquation → DNA guanine N7 binding → intrastrand cross-link",
            enzyme_families=["(not enzymatic — cisplatin/carboplatin pharmaceutical)"],
            mineral_analogs=["Sperrylite (PtAs₂)", "Cooperite (PtS)"],
            eh_window=(0.2, 1.2),
            ph_window=(3.0, 8.0),
            notes="No known Pt enzymes. Pt complexes are pharmaceutical (cisplatin) "
                  "or industrial catalysts. Square planar d⁸ is kinetically labile for ligand exchange.",
        ),
    },
    # ── Lanthanides ──
    {
        "elements": ["La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy",
                      "Ho", "Er", "Tm", "Yb", "Lu", "Y"],
        "ox_states": [3],
        "geometries": None,  # any geometry
        "donor_pattern": lambda ds: True,
        "prediction": CatalyticPrediction(
            function="Methanol oxidation / Lewis acid catalysis (lanthanome)",
            confidence="moderate",
            mechanism="LnIII activates methanol dehydrogenase (XoxF-type MDH); "
                      "Lewis acid polarization of substrate",
            enzyme_families=["XoxF-type methanol dehydrogenase (MDH)", "Lanmodulin",
                             "ExaF ethanol dehydrogenase"],
            mineral_analogs=["Monazite ((Ce,La)PO₄)", "Bastnäsite ((Ce,La)(CO₃)F)",
                             "Xenotime (YPO₄)"],
            eh_window=(-0.2, 0.5),
            ph_window=(5.0, 8.0),
            notes="Lanthanide biochemistry discovered ~2011. Ln³⁺ replaces Ca²⁺ in MDH. "
                  "Higher charge density → stronger Lewis acid → better methanol activation. "
                  "Lanmodulin (LanM) selectively binds Ln over Ca by 10⁸-fold.",
        ),
    },
    # ── Tungsten ──
    {
        "elements": ["W"],
        "ox_states": [4, 6],
        "geometries": ["square pyramidal", "distorted square pyramidal", "octahedral"],
        "donor_pattern": lambda ds: "S" in ds or "O" in ds,
        "prediction": CatalyticPrediction(
            function="Oxo-transfer (anaerobic, low-potential)",
            confidence="moderate",
            mechanism="Analogous to Mo oxo-transfer but at lower reduction potential. "
                      "WVI=O + substrate → WIV + substrate-O",
            enzyme_families=["Aldehyde ferredoxin oxidoreductase", "Formaldehyde ferredoxin oxidoreductase",
                             "W-containing formate dehydrogenase"],
            mineral_analogs=["Wolframite ((Fe,Mn)WO₄)", "Scheelite (CaWO₄)", "Tungstenite (WS₂)"],
            eh_window=(-0.5, 0.2),
            ph_window=(5.0, 9.0),
            notes="W enzymes are restricted to anaerobic archaea/bacteria. "
                  "Lower reduction potential than Mo analogs → operates in more reducing environments.",
        ),
    },
]


# ─────────────────────────────────────────────────────────────────────
# Prediction engine
# ─────────────────────────────────────────────────────────────────────

def predict_function(
    site: CoordinationSite,
    oxidation_state: int = 2,
) -> list[CatalyticPrediction]:
    """
    Predict catalytic function(s) from coordination site properties.

    Returns a ranked list of predictions (best match first).
    """
    element = site.metal.element
    geo = (site.geometry or "").lower()
    donor_set = site.donor_set

    matches: list[tuple[int, CatalyticPrediction]] = []

    for rule in _RULES:
        score = 0

        # Element match (required)
        if element not in rule["elements"]:
            continue
        score += 10

        # Oxidation state match
        if oxidation_state in rule["ox_states"]:
            score += 5

        # Geometry match
        if rule["geometries"] is None:
            score += 3  # any geometry
        else:
            for g in rule["geometries"]:
                if g in geo or geo in g:
                    score += 5
                    break

        # Donor set match
        if rule["donor_pattern"](donor_set):
            score += 3

        matches.append((score, rule["prediction"]))

    # Sort by score descending
    matches.sort(key=lambda x: -x[0])
    return [pred for _, pred in matches]


def entatic_state_analysis(site: CoordinationSite) -> Optional[dict]:
    """
    Analyze whether the site shows entatic state characteristics.

    The entatic state hypothesis: proteins distort metal coordination
    away from the thermodynamic preference of the free ion to optimize
    catalytic function (usually electron transfer rate).

    Returns a dict with:
        - preferred_geometry: what the free ion prefers
        - actual_geometry: what the protein enforces
        - distortion_type: description of the strain
        - entatic_score: 0-1 (0 = relaxed, 1 = highly strained)
        - functional_benefit: why the protein does this
    """
    from terrahedral.analysis import tau4, tau5

    element = site.metal.element
    geo = (site.geometry or "").lower()
    cn = site.coordination_number

    # Free ion geometry preferences (thermodynamic minimum)
    ION_PREFERENCES = {
        ("Cu", 2): {"preferred": "square planar", "reason": "d⁹ Jahn-Teller, tetragonal elongation"},
        ("Cu", 1): {"preferred": "linear or trigonal", "reason": "d¹⁰ no LFSE preference, small ion"},
        ("Fe", 2): {"preferred": "octahedral", "reason": "d⁶ maximizes LFSE in Oh"},
        ("Fe", 3): {"preferred": "octahedral", "reason": "d⁵ HS has zero LFSE but octahedral by size"},
        ("Zn", 2): {"preferred": "tetrahedral", "reason": "d¹⁰ no LFSE, small size favors Td"},
        ("Ni", 2): {"preferred": "square planar or octahedral", "reason": "d⁸ strong LFSE in both"},
        ("Co", 2): {"preferred": "tetrahedral or octahedral", "reason": "d⁷ LFSE modest, geometry flexible"},
        ("Co", 3): {"preferred": "octahedral", "reason": "d⁶ LS maximizes LFSE in Oh"},
        ("Mn", 2): {"preferred": "octahedral", "reason": "d⁵ HS zero LFSE, octahedral by size"},
        ("Pt", 2): {"preferred": "square planar", "reason": "d⁸ strong-field 3rd row, always sq. planar"},
        ("Pd", 2): {"preferred": "square planar", "reason": "d⁸ strong-field 2nd row"},
    }

    result = {
        "preferred_geometry": None,
        "actual_geometry": geo,
        "distortion_type": None,
        "entatic_score": 0.0,
        "functional_benefit": None,
    }

    # Try multiple oxidation states
    for ox in [2, 3, 1, 4]:
        key = (element, ox)
        if key in ION_PREFERENCES:
            pref = ION_PREFERENCES[key]
            result["preferred_geometry"] = pref["preferred"]

            # Calculate entatic score based on geometry mismatch
            preferred = pref["preferred"].lower()
            actual = geo

            # Direct match → no strain
            if preferred in actual or actual in preferred:
                result["entatic_score"] = 0.1
                result["distortion_type"] = "Minor distortion from ideal"
                result["functional_benefit"] = "Near-optimal geometry; minimal entatic effect"
                break

            # Classic entatic case: Cu²⁺ in tetrahedral (should be sq. planar)
            if element == "Cu" and ox == 2 and "tetra" in actual:
                t4_val = tau4(site)
                result["entatic_score"] = 0.9
                result["distortion_type"] = (
                    f"CuII forced into tetrahedral (prefers square planar). "
                    f"τ₄ = {t4_val:.3f}" if t4_val else "CuII forced into tetrahedral"
                )
                result["functional_benefit"] = (
                    "Minimizes reorganization energy for CuII/CuI electron transfer. "
                    "Tetrahedral is ideal for CuI (d¹⁰) — protein pre-organizes the "
                    "oxidized state geometry to match the reduced state, accelerating ET by ~10³."
                )
                break

            # Geometry swap
            if "tetra" in preferred and "planar" in actual:
                result["entatic_score"] = 0.7
                result["distortion_type"] = f"Ion prefers tetrahedral but protein enforces {actual}"
                result["functional_benefit"] = "Geometry strain may tune reduction potential or substrate access"
            elif "planar" in preferred and "tetra" in actual:
                result["entatic_score"] = 0.8
                result["distortion_type"] = f"Ion prefers square planar but protein enforces {actual}"
                result["functional_benefit"] = "May optimize electron transfer rate by reducing reorganization energy"
            elif "octa" in preferred and ("tetra" in actual or "pyram" in actual):
                result["entatic_score"] = 0.6
                result["distortion_type"] = f"Ion prefers octahedral but site is {actual} (CN={cn})"
                result["functional_benefit"] = "Reduced CN may open coordination site for substrate binding"
            else:
                result["entatic_score"] = 0.4
                result["distortion_type"] = f"Moderate mismatch: prefers {preferred}, actual {actual}"
                result["functional_benefit"] = "Geometry distortion may tune reactivity"
            break

    return result if result["preferred_geometry"] else None

"""
Eh-pH stability windows for metal species.

Maps metal + oxidation state → Eh-pH boundary conditions for
aqueous stability.  Used to compare metalloenzyme operating
conditions with mineral stability fields.

Eh = reduction potential (V vs SHE)
pH = -log[H+]

Regions are defined as *polygons* in (pH, Eh) space so that
Nernst-equation diagonal boundaries render correctly.  Legacy
rectangular regions (eh_min/max, ph_min/max) are still accepted;
the renderer promotes them to four-corner polygons automatically.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


@dataclass
class EhPhRegion:
    """An Eh-pH stability region for a metal species.

    Define either:
      * vertices – list of (pH, Eh) polygon corners (preferred), or
      * eh_min/max + ph_min/max – axis-aligned rectangle (legacy).
    If vertices is non-empty it takes precedence.
    """
    species: str            # e.g. "Fe2+(aq)", "FeOOH(s)"
    phase: str              # "aqueous", "solid", "gaseous"
    eh_min: float = 0.0     # V vs SHE  (used only if vertices empty)
    eh_max: float = 0.0
    ph_min: float = 0.0
    ph_max: float = 0.0
    color: str = "#888"     # for rendering
    label: str = ""         # display label
    notes: str = ""
    vertices: list[tuple[float, float]] = field(default_factory=list)
    # vertices are (pH, Eh) pairs

    def polygon(self) -> list[tuple[float, float]]:
        """Return polygon vertices, promoting rect to polygon if needed."""
        if self.vertices:
            return list(self.vertices)
        return [
            (self.ph_min, self.eh_min),
            (self.ph_max, self.eh_min),
            (self.ph_max, self.eh_max),
            (self.ph_min, self.eh_max),
        ]


@dataclass
class EhPhDiagram:
    """Complete Eh-pH diagram for an element."""
    element: str
    regions: list[EhPhRegion] = field(default_factory=list)
    water_oxidation_line: bool = True  # show O2/H2O boundary
    water_reduction_line: bool = True  # show H2/H2O boundary
    biological_window: bool = True     # show typical biological range


# -----------------------------------------------------------------
# Helper: Nernst boundary
# Eh = E0 - (0.0592 * m / n) * pH   where m = H+, n = e-
# -----------------------------------------------------------------

def _nernst(e0: float, m: int, n: int, ph: float) -> float:
    return e0 - 0.0592 * m / n * ph


# =================================================================
#  Eh-pH data for key metals
#  Polygonal stability fields approximated from published Pourbaix
#  diagrams (25 deg C, 10^-6 M dissolved species, 1 atm gas).
# =================================================================

DIAGRAMS: dict[str, EhPhDiagram] = {}


# -- Vanadium -------------------------------------------------------
# Key couples:
#   V3+ + e- -> V2+                E0 = -0.255 V
#   VO2+ + 2H+ + e- -> V3+ + H2O  E0 ~  0.337   slope -0.118
#   VO2+ + 2H+ + e- -> VO2+ + H2O E0 ~  1.00    slope -0.118
#   V2O5 dissolves -> VO2+ (acid) or HVO42- (alkaline)

DIAGRAMS["V"] = EhPhDiagram(
    element="V",
    regions=[
        # V2+(aq) -- very reducing, acidic (pH 0-4)
        EhPhRegion(
            "V2+(aq)", "aqueous",
            color="#5A8FD4", label="V\u00b2\u207a",
            notes="Vanadous ion. Violet. Strong reductant.",
            vertices=[
                (0.0, -1.0), (4.0, -1.0), (4.0, -0.255), (0.0, -0.255),
            ],
        ),
        # V3+(aq) -- reducing, acidic (pH 0-4)
        EhPhRegion(
            "V3+(aq)", "aqueous",
            color="#1D9E75", label="V\u00b3\u207a",
            notes="Vanadic ion. Green. Moderately reducing.",
            vertices=[
                (0.0, -0.255), (4.0, -0.255),
                (4.0, -0.1366), (0.0, 0.337),
            ],
        ),
        # VO2+(aq) -- vanadyl, mildly oxidizing, acidic (pH 0-4)
        EhPhRegion(
            "VO2+(aq)", "aqueous",
            color="#D4A85A", label="VO\u00b2\u207a",
            notes="Vanadyl ion. Blue. Dominant V(IV) in acidic solution.",
            vertices=[
                (0.0, 0.337), (4.0, -0.1366),
                (4.0, 0.5264), (0.0, 1.0),
            ],
        ),
        # VO2+(aq) -- pervanadyl, strongly oxidizing, acidic (pH 0-4)
        EhPhRegion(
            "VO2+(aq)", "aqueous",
            color="#D85A30", label="VO\u2082\u207a",
            notes="Pervanadyl. Yellow. V(V) in strongly acidic solution.",
            vertices=[
                (0.0, 1.0), (4.0, 0.5264),
                (4.0, 1.4), (0.0, 1.4),
            ],
        ),
        # V2O3(s) -- sesquioxide (pH 4-9, reducing)
        EhPhRegion(
            "V2O3(s)", "solid",
            color="#7A8B99", label="V\u2082O\u2083",
            notes="Karelianite. Black. Stable under reducing conditions.",
            vertices=[
                (4.0, -1.0), (9.0, -1.0),
                (9.0, -0.7286000000000001), (4.0, -0.1366),
            ],
        ),
        # VO2(s) -- paramontroseite (pH 4-9, intermediate)
        EhPhRegion(
            "VO2(s)", "solid",
            color="#534AB7", label="VO\u2082",
            notes="Paramontroseite. V(IV) oxide.",
            vertices=[
                (4.0, -0.1366), (9.0, -0.7286000000000001),
                (9.0, -0.0656000000000001), (4.0, 0.5264),
            ],
        ),
        # V2O5(s) -- pentoxide (pH 4-9, oxidizing)
        EhPhRegion(
            "V2O5(s)", "solid",
            color="#D4537E", label="V\u2082O\u2085",
            notes="Shcherbinaite. Orange. Industrial oxidation catalyst.",
            vertices=[
                (4.0, 0.5264), (9.0, -0.0656000000000001),
                (9.0, 1.4), (4.0, 1.4),
            ],
        ),
        # V(OH)3(s) -- hydroxide, alkaline reducing (pH 9-14)
        EhPhRegion(
            "V(OH)3(s)", "solid",
            color="#6AAAB9", label="V(OH)\u2083",
            notes="V(III) hydroxide. Precipitates at alkaline pH.",
            vertices=[
                (9.0, -1.0), (14.0, -1.0),
                (14.0, -0.6576), (9.0, -0.0656000000000001),
                (9.0, -0.7286000000000001),
            ],
        ),
        # HVO42-(aq) -- vanadate, alkaline oxidizing (pH 9-14)
        EhPhRegion(
            "HVO42-(aq)", "aqueous",
            color="#8BB8C8", label="HVO\u2084\u00b2\u207b",
            notes="Vanadate. Dominant V(V) at alkaline pH. Bioavailable.",
            vertices=[
                (9.0, -0.0656000000000001), (14.0, -0.6576),
                (14.0, 1.4), (9.0, 1.4),
            ],
        ),
    ],
)


# -- Iron -----------------------------------------------------------
DIAGRAMS["Fe"] = EhPhDiagram(
    element="Fe",
    regions=[
        EhPhRegion("Fe2+(aq)", "aqueous",
                   color="#5A8FD4", label="Fe\u00b2\u207a",
                   notes="Soluble ferrous iron. Mobile in reducing, acidic conditions.",
                   vertices=[
                       (0.0, -0.62), (6.5, -0.62),
                       (6.5, _nernst(0.771, 0, 1, 0.0)),
                       (0.0, _nernst(0.771, 0, 1, 0.0)),
                   ]),
        EhPhRegion("Fe3+(aq)", "aqueous",
                   color="#D85A30", label="Fe\u00b3\u207a",
                   notes="Soluble ferric. Only at very low pH.",
                   vertices=[
                       (0.0, _nernst(0.771, 0, 1, 0.0)),
                       (3.0, _nernst(0.771, 0, 1, 0.0)),
                       (3.0, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Fe(OH)3 / FeOOH", "solid",
                   color="#D4A85A", label="FeOOH(s)",
                   notes="Goethite/ferrihydrite. Insoluble above pH ~3.",
                   vertices=[
                       (3.0, _nernst(0.771, 0, 1, 0.0)),
                       (3.0, 1.4), (14.0, 1.4),
                       (14.0, _nernst(-0.044, 1, 1, 14.0)),
                       (6.5, _nernst(-0.044, 1, 1, 6.5)),
                   ]),
        EhPhRegion("Fe(OH)2", "solid",
                   color="#7A8B99", label="Fe(OH)\u2082(s)",
                   notes="Ferrous hydroxide. Metastable.",
                   vertices=[
                       (6.5, -0.62), (14.0, -0.62),
                       (14.0, _nernst(-0.044, 1, 1, 14.0)),
                       (6.5, _nernst(-0.044, 1, 1, 6.5)),
                   ]),
        EhPhRegion("Fe(s)", "solid",
                   color="#2A2D35", label="Fe(s)",
                   notes="Metallic iron. Very reducing.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, -0.62), (0.0, -0.62),
                   ]),
    ],
)


# -- Copper ---------------------------------------------------------
DIAGRAMS["Cu"] = EhPhDiagram(
    element="Cu",
    regions=[
        EhPhRegion("Cu2+(aq)", "aqueous",
                   color="#378ADD", label="Cu\u00b2\u207a",
                   notes="Soluble cupric ion.",
                   vertices=[
                       (0.0, 0.34), (5.5, 0.34),
                       (5.5, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("CuO", "solid",
                   color="#1D9E75", label="CuO(s)",
                   notes="Tenorite. Stable cupric oxide.",
                   vertices=[
                       (5.5, 0.34), (5.5, 1.4), (14.0, 1.4),
                       (14.0, _nernst(0.1, 2, 2, 14.0)),
                   ]),
        EhPhRegion("Cu2O", "solid",
                   color="#D85A30", label="Cu\u2082O(s)",
                   notes="Cuprite. Low-Eh copper oxide.",
                   vertices=[
                       (5.5, 0.34),
                       (14.0, _nernst(0.1, 2, 2, 14.0)),
                       (14.0, _nernst(-0.36, 2, 2, 14.0)),
                       (7.0, -0.12),
                   ]),
        EhPhRegion("Cu(s)", "solid",
                   color="#888780", label="Cu(s)",
                   notes="Native copper. Very reducing.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, _nernst(-0.36, 2, 2, 14.0)),
                       (7.0, -0.12),
                       (5.5, 0.34), (0.0, 0.34),
                   ]),
    ],
)


# -- Manganese ------------------------------------------------------
DIAGRAMS["Mn"] = EhPhDiagram(
    element="Mn",
    regions=[
        EhPhRegion("Mn2+(aq)", "aqueous",
                   color="#D4537E", label="Mn\u00b2\u207a",
                   notes="Soluble. Dominant species in reducing conditions.",
                   vertices=[
                       (0.0, -1.0), (9.0, -1.0),
                       (9.0, _nernst(0.0, 2, 2, 9.0)),
                       (0.0, 0.6),
                   ]),
        EhPhRegion("MnO2", "solid",
                   color="#993556", label="MnO\u2082(s)",
                   notes="Pyrolusite/birnessite. Natural water oxidation catalyst.",
                   vertices=[
                       (0.0, 0.6),
                       (0.0, 1.4), (14.0, 1.4),
                       (14.0, _nernst(1.23, 4, 2, 14.0)),
                       (9.0, _nernst(0.0, 2, 2, 9.0)),
                   ]),
        EhPhRegion("Mn(OH)2", "solid",
                   color="#7A8B99", label="Mn(OH)\u2082(s)",
                   notes="Pyrochroite. Reduced Mn hydroxide.",
                   vertices=[
                       (9.0, -1.0), (14.0, -1.0),
                       (14.0, _nernst(1.23, 4, 2, 14.0)),
                       (9.0, _nernst(0.0, 2, 2, 9.0)),
                   ]),
    ],
)


# -- Zinc -----------------------------------------------------------
DIAGRAMS["Zn"] = EhPhDiagram(
    element="Zn",
    regions=[
        EhPhRegion("Zn2+(aq)", "aqueous",
                   color="#7F77DD", label="Zn\u00b2\u207a",
                   notes="Soluble. Wide stability range.",
                   vertices=[
                       (0.0, -1.0), (8.5, -1.0),
                       (8.5, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Zn(OH)2 / ZnO", "solid",
                   color="#534AB7", label="ZnO(s)",
                   notes="Zincite. Amphoteric.",
                   vertices=[
                       (8.5, -1.0), (12.5, -1.0),
                       (12.5, 1.4), (8.5, 1.4),
                   ]),
        EhPhRegion("ZnO22-(aq)", "aqueous",
                   color="#B4B2A9", label="ZnO\u2082\u00b2\u207b",
                   notes="Zincate ion. High pH.",
                   vertices=[
                       (12.5, -1.0), (14.0, -1.0),
                       (14.0, 1.4), (12.5, 1.4),
                   ]),
    ],
)


# -- Molybdenum -----------------------------------------------------
DIAGRAMS["Mo"] = EhPhDiagram(
    element="Mo",
    regions=[
        EhPhRegion("Mo(s)", "solid",
                   color="#2A2D35", label="Mo(s)",
                   notes="Metallic Mo. Very reducing.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, -0.9007999999999999), (0.0, -0.072),
                   ]),
        EhPhRegion("MoO2(s)", "solid",
                   color="#7A8B99", label="MoO\u2082(s)",
                   notes="Reduced Mo oxide.",
                   vertices=[
                       (0.0, -0.072), (14.0, -0.9007999999999999),
                       (14.0, -0.8288), (0.0, -0.0),
                   ]),
        EhPhRegion("MoO42-(aq)", "aqueous",
                   color="#7F77DD", label="MoO\u2084\u00b2\u207b",
                   notes="Molybdate. Bioavailable form.",
                   vertices=[
                       (0.0, -0.0), (14.0, -0.8288),
                       (14.0, 1.4), (0.0, 1.4),
                   ]),
    ],
)


# -- Nickel ---------------------------------------------------------
DIAGRAMS["Ni"] = EhPhDiagram(
    element="Ni",
    regions=[
        EhPhRegion("Ni2+(aq)", "aqueous",
                   color="#97C459", label="Ni\u00b2\u207a",
                   notes="Soluble nickel(II).",
                   vertices=[
                       (0.0, _nernst(-0.257, 0, 2, 0.0)),
                       (8.0, _nernst(-0.257, 0, 2, 0.0)),
                       (8.0, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Ni(OH)2", "solid",
                   color="#3B6D11", label="Ni(OH)\u2082(s)",
                   notes="Theophrastite.",
                   vertices=[
                       (8.0, _nernst(-0.257, 0, 2, 0.0)),
                       (14.0, _nernst(-0.257, 0, 2, 0.0)),
                       (14.0, _nernst(0.49, 2, 1, 14.0)),
                       (8.0, _nernst(0.49, 2, 1, 8.0)),
                   ]),
        EhPhRegion("NiOOH", "solid",
                   color="#D4A85A", label="NiOOH(s)",
                   notes="Nickel oxyhydroxide.",
                   vertices=[
                       (8.0, _nernst(0.49, 2, 1, 8.0)),
                       (14.0, _nernst(0.49, 2, 1, 14.0)),
                       (14.0, 1.4), (8.0, 1.4),
                   ]),
        EhPhRegion("Ni(s)", "solid",
                   color="#5F5E5A", label="Ni(s)",
                   notes="Metallic nickel.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, _nernst(-0.257, 0, 2, 0.0)),
                       (0.0, _nernst(-0.257, 0, 2, 0.0)),
                   ]),
    ],
)


# -- Cobalt ---------------------------------------------------------
DIAGRAMS["Co"] = EhPhDiagram(
    element="Co",
    regions=[
        EhPhRegion("Co2+(aq)", "aqueous",
                   color="#1D9E75", label="Co\u00b2\u207a",
                   notes="Soluble cobalt(II). Rose-pink.",
                   vertices=[
                       (0.0, _nernst(-0.277, 0, 2, 0.0)),
                       (7.5, _nernst(-0.277, 0, 2, 0.0)),
                       (7.5, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Co(OH)2", "solid",
                   color="#0F6E56", label="Co(OH)\u2082(s)",
                   notes="Cobalt hydroxide.",
                   vertices=[
                       (7.5, _nernst(-0.277, 0, 2, 0.0)),
                       (14.0, _nernst(-0.277, 0, 2, 0.0)),
                       (14.0, _nernst(0.17, 1, 1, 14.0)),
                       (7.5, _nernst(0.17, 1, 1, 7.5)),
                   ]),
        EhPhRegion("Co3O4 / CoOOH", "solid",
                   color="#D85A30", label="Co\u2083O\u2084(s)",
                   notes="Cobalt spinel. Mixed-valence.",
                   vertices=[
                       (7.5, _nernst(0.17, 1, 1, 7.5)),
                       (14.0, _nernst(0.17, 1, 1, 14.0)),
                       (14.0, 1.4), (7.5, 1.4),
                   ]),
        EhPhRegion("Co(s)", "solid",
                   color="#2A2D35", label="Co(s)",
                   notes="Metallic cobalt.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, _nernst(-0.277, 0, 2, 0.0)),
                       (0.0, _nernst(-0.277, 0, 2, 0.0)),
                   ]),
    ],
)


# -- Tungsten -------------------------------------------------------
DIAGRAMS["W"] = EhPhDiagram(
    element="W",
    regions=[
        EhPhRegion("W(s)", "solid",
                   color="#2A2D35", label="W(s)",
                   notes="Metallic tungsten. Very reducing.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, -0.9188), (0.0, -0.09),
                   ]),
        EhPhRegion("WO2(s)", "solid",
                   color="#5F5E5A", label="WO\u2082(s)",
                   notes="Brown oxide.",
                   vertices=[
                       (0.0, -0.09), (14.0, -0.9188),
                       (14.0, -0.7787999999999999), (0.0, 0.05),
                   ]),
        EhPhRegion("WO3(s)", "solid",
                   color="#BA7517", label="WO\u2083(s)",
                   notes="Tungstic oxide. Yellow.",
                   vertices=[
                       (0.0, 0.05), (6.0, -0.3052),
                       (6.0, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("WO42-(aq)", "aqueous",
                   color="#888780", label="WO\u2084\u00b2\u207b",
                   notes="Tungstate. Bioavailable. Used by tungstoenzymes.",
                   vertices=[
                       (6.0, -0.3052), (14.0, -0.7787999999999999),
                       (14.0, 1.4), (6.0, 1.4),
                   ]),
    ],
)


# -- Chromium -------------------------------------------------------
DIAGRAMS["Cr"] = EhPhDiagram(
    element="Cr",
    regions=[
        EhPhRegion("Cr3+(aq)", "aqueous",
                   color="#BA7517", label="Cr\u00b3\u207a",
                   notes="Chromic ion. Green.",
                   vertices=[
                       (0.0, _nernst(-0.74, 0, 1, 0.0)),
                       (4.0, _nernst(-0.74, 0, 1, 0.0)),
                       (4.0, 0.6), (0.0, 0.6),
                   ]),
        EhPhRegion("Cr2O3 / Cr(OH)3", "solid",
                   color="#7A8B99", label="Cr\u2082O\u2083(s)",
                   notes="Eskolaite/chromia. Very stable. Passivation layer.",
                   vertices=[
                       (4.0, -0.75), (12.0, -0.75),
                       (12.0, _nernst(0.6, 2, 1, 12.0)),
                       (4.0, _nernst(0.6, 2, 1, 4.0)),
                   ]),
        EhPhRegion("CrO42-(aq)", "aqueous",
                   color="#D85A30", label="CrO\u2084\u00b2\u207b",
                   notes="Chromate. Cr(VI). Toxic, strong oxidant.",
                   vertices=[
                       (0.0, 0.6), (4.0, 0.6),
                       (4.0, _nernst(0.6, 2, 1, 4.0)),
                       (12.0, _nernst(0.6, 2, 1, 12.0)),
                       (14.0, _nernst(0.6, 2, 1, 14.0)),
                       (14.0, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Cr(s)", "solid",
                   color="#2A2D35", label="Cr(s)",
                   notes="Metallic chromium.",
                   vertices=[
                       (0.0, -1.0), (14.0, -1.0),
                       (14.0, -0.75), (12.0, -0.75),
                       (4.0, -0.75),
                       (4.0, _nernst(-0.74, 0, 1, 0.0)),
                       (0.0, _nernst(-0.74, 0, 1, 0.0)),
                   ]),
    ],
)


# -- Neodymium (representative lanthanide) --------------------------
DIAGRAMS["Nd"] = EhPhDiagram(
    element="Nd",
    regions=[
        EhPhRegion("Nd3+(aq)", "aqueous",
                   color="#8BB8C8", label="Nd\u00b3\u207a",
                   notes="Soluble Nd. Only +3 in aqueous.",
                   vertices=[
                       (0.0, -1.0), (7.5, -1.0),
                       (7.5, 1.4), (0.0, 1.4),
                   ]),
        EhPhRegion("Nd(OH)3", "solid",
                   color="#5A9AA9", label="Nd(OH)\u2083(s)",
                   notes="Insoluble Ln hydroxide.",
                   vertices=[
                       (7.5, -1.0), (14.0, -1.0),
                       (14.0, 1.4), (7.5, 1.4),
                   ]),
    ],
)

# Apply to all lanthanides
for _el in ["La", "Ce", "Pr", "Sm", "Eu", "Gd", "Tb", "Dy",
            "Ho", "Er", "Tm", "Yb", "Lu", "Y"]:
    DIAGRAMS[_el] = EhPhDiagram(
        element=_el,
        regions=[
            EhPhRegion(f"{_el}3+(aq)", "aqueous",
                       color="#8BB8C8", label=f"{_el}\u00b3\u207a",
                       vertices=[
                           (0.0, -1.0), (7.5, -1.0),
                           (7.5, 1.4), (0.0, 1.4),
                       ]),
            EhPhRegion(f"{_el}(OH)3", "solid",
                       color="#5A9AA9", label=f"{_el}(OH)\u2083(s)",
                       vertices=[
                           (7.5, -1.0), (14.0, -1.0),
                           (14.0, 1.4), (7.5, 1.4),
                       ]),
        ],
    )

# Ce: extra +4 state
DIAGRAMS["Ce"].regions.insert(1, EhPhRegion(
    "CeO2", "solid",
    color="#3A7A99", label="CeO\u2082(s)",
    notes="Ceria. CeIII->CeIV. Only Ln with accessible +4 in aqueous.",
    vertices=[
        (2.0, _nernst(1.72, 4, 1, 2.0)),
        (14.0, _nernst(1.72, 4, 1, 14.0)),
        (14.0, 1.4), (2.0, 1.4),
    ],
))

# Eu: +2 state
DIAGRAMS["Eu"].regions.insert(0, EhPhRegion(
    "Eu2+(aq)", "aqueous",
    color="#6AAAB9", label="Eu\u00b2\u207a",
    notes="Divalent Eu. Reducing conditions. Eu anomaly.",
    vertices=[
        (0.0, -1.0), (7.0, -1.0),
        (7.0, -0.35), (0.0, -0.35),
    ],
))


# -----------------------------------------------------------------
# Water stability lines (for rendering)
# -----------------------------------------------------------------

def water_oxidation_eh(ph: float) -> float:
    """O2/H2O line: Eh = 1.229 - 0.0592 * pH  (25 C, 1 atm O2)."""
    return 1.229 - 0.0592 * ph


def water_reduction_eh(ph: float) -> float:
    """H2/H2O line: Eh = -0.0592 * pH  (25 C, 1 atm H2)."""
    return -0.0592 * ph


# -----------------------------------------------------------------
# Biological windows
# -----------------------------------------------------------------

BIOLOGICAL_WINDOWS = {
    "cytoplasm": {"eh_range": (-0.4, -0.2), "ph_range": (6.8, 7.4),
                  "label": "Cytoplasm", "color": "rgba(90,196,126,0.15)"},
    "mitochondria": {"eh_range": (-0.3, 0.0), "ph_range": (7.0, 8.0),
                     "label": "Mitochondria", "color": "rgba(90,143,212,0.15)"},
    "periplasm": {"eh_range": (-0.2, 0.2), "ph_range": (6.5, 7.5),
                  "label": "Periplasm", "color": "rgba(212,168,90,0.15)"},
    "acidic_vent": {"eh_range": (-0.4, 0.0), "ph_range": (2.0, 5.0),
                    "label": "Acidic hydrothermal", "color": "rgba(212,90,90,0.15)"},
    "alkaline_vent": {"eh_range": (-0.6, -0.2), "ph_range": (9.0, 11.0),
                      "label": "Alkaline vent", "color": "rgba(143,90,212,0.15)"},
    "soil_oxic": {"eh_range": (0.2, 0.8), "ph_range": (4.0, 8.0),
                  "label": "Oxic soil", "color": "rgba(196,180,90,0.15)"},
    "soil_anoxic": {"eh_range": (-0.3, 0.1), "ph_range": (5.0, 8.0),
                    "label": "Anoxic soil/sediment", "color": "rgba(122,139,153,0.15)"},
    "seawater": {"eh_range": (0.3, 0.5), "ph_range": (7.8, 8.3),
                 "label": "Seawater", "color": "rgba(90,143,212,0.15)"},
}


def get_diagram(element: str) -> Optional[EhPhDiagram]:
    """Return the Eh-pH diagram for an element, or None."""
    return DIAGRAMS.get(element)


def get_operating_window(site, oxidation_state: int = 2) -> Optional[dict]:
    """
    Estimate the Eh-pH operating window for a coordination site based on
    its metal, oxidation state, and predicted catalytic function.

    Returns a dict with eh_range, ph_range, function, and examples,
    or None if no prediction can be made.
    """
    from terrahedral.analysis.catalysis import infer_function

    predictions = infer_function(site)
    if not predictions:
        return None
    top = predictions[0]

    # Estimate Eh-pH window from the metal's Pourbaix diagram
    el = site.metal.element
    diagram = get_diagram(el)

    # Default biological window
    eh_range = (-0.4, 0.4)
    ph_range = (5.0, 9.0)

    if diagram and diagram.regions:
        # Find the region(s) matching the oxidation state
        # Use a reasonable biological range as default
        for region in diagram.regions:
            if region.phase == "aqueous" and f"{oxidation_state}+" in region.species:
                verts = region.polygon()
                if verts:
                    phs = [v[0] for v in verts]
                    ehs = [v[1] for v in verts]
                    eh_range = (min(ehs), max(ehs))
                    ph_range = (min(phs), max(phs))
                break

    return {
        "eh_range": eh_range,
        "ph_range": ph_range,
        "function": top.function,
        "examples": top.examples,
    }

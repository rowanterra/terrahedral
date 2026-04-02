"""
terrahedral — CIF -> coordination geometry analysis
====================================================

Extract metal coordination environments from CIF files and
analyze geometry, bonding, and catalytic properties.

    >>> from terrahedral import CoordinationSite
    >>> site = CoordinationSite.from_cif("structure.cif", metal="Cu1")
    >>> site.classify_geometry()
    >>> print(site.summary())

Supports both small-molecule CIF (COD/CCDC) and macromolecular
mmCIF (PDB) formats.
"""

__version__ = "0.1.0"

from terrahedral.core import CoordinationSite, Atom, Bond, Angle, METAL_ELEMENTS
from terrahedral.io import load_cif, load_mmcif
from terrahedral.fetch import fetch_pdb, fetch_cod, search_pdb, search_cod
from terrahedral.catalytic import CatalyticFrame, CatalyticPathway
from terrahedral.parsers.xyz_parser import load_xyz
from terrahedral.analysis.catalysis import (
    CatalyticCycle, CycleStep,
    infer_function, analyze_entatic,
    FunctionPrediction, EntaicAnalysis,
)
from terrahedral.analysis.alignment import compare_shells, ShellComparison

__all__ = [
    "CoordinationSite",
    "Atom",
    "Bond",
    "Angle",
    "METAL_ELEMENTS",
    "load_cif",
    "load_mmcif",
    "load_xyz",
    "fetch_pdb",
    "fetch_cod",
    "search_pdb",
    "search_cod",
    "CatalyticFrame",
    "CatalyticPathway",
    "CatalyticCycle",
    "CycleStep",
    "infer_function",
    "analyze_entatic",
    "FunctionPrediction",
    "EntaicAnalysis",
    "compare_shells",
    "ShellComparison",
]

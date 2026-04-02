"""Tests for terrahedral (barebones — no rendering)."""

import math
from pathlib import Path

import pytest

from terrahedral.core import Atom, Bond, Angle, CoordinationSite
from terrahedral.transforms import frac_to_cart, distance, angle_at_center, rotate_y
from terrahedral.analysis import tau4, tau5, classify


# ── Transform tests ─────────────────────────────────────────────────────────

def test_frac_to_cart_cubic():
    """Cubic cell: fractional (0.5, 0.5, 0.5) -> (5, 5, 5) for a=10."""
    x, y, z = frac_to_cart(0.5, 0.5, 0.5, 10, 10, 10, 90, 90, 90)
    assert abs(x - 5.0) < 1e-6
    assert abs(y - 5.0) < 1e-6
    assert abs(z - 5.0) < 1e-6


def test_frac_to_cart_monoclinic():
    """Monoclinic cell: beta != 90 should shift x."""
    x1, y1, z1 = frac_to_cart(0, 0, 1.0, 10, 10, 10, 90, 90, 90)
    x2, y2, z2 = frac_to_cart(0, 0, 1.0, 10, 10, 10, 90, 120, 90)
    # In monoclinic, c-axis has an x-component
    assert abs(x1 - 0.0) < 1e-6  # orthorhombic: no x shift
    assert x2 != 0.0              # monoclinic: x shifts


def test_distance():
    assert abs(distance((0, 0, 0), (3, 4, 0)) - 5.0) < 1e-10


def test_angle_at_center():
    # 90 degree angle
    a = angle_at_center((1, 0, 0), (0, 0, 0), (0, 1, 0))
    assert abs(a - 90.0) < 1e-6

    # 180 degree angle
    a = angle_at_center((1, 0, 0), (0, 0, 0), (-1, 0, 0))
    assert abs(a - 180.0) < 1e-6


def test_rotate_y():
    x, y, z = rotate_y((1, 0, 0), 90)
    assert abs(x) < 1e-6
    assert abs(z - (-1.0)) < 1e-6


# ── Core data model tests ──────────────────────────────────────────────────

def test_atom_distance():
    a1 = Atom("Cu", "Cu", 0, 0, 0)
    a2 = Atom("N1", "N", 2, 0, 0)
    assert abs(a1.distance_to(a2) - 2.0) < 1e-10


def test_atom_centered():
    origin = Atom("Cu", "Cu", 5, 5, 5)
    a = Atom("N1", "N", 7, 5, 5)
    c = a.centered_on(origin)
    assert abs(c.x - 2.0) < 1e-10
    assert abs(c.y) < 1e-10


def test_coordination_site_cn():
    site = CoordinationSite(
        metal=Atom("Cu", "Cu", 0, 0, 0),
        ligands=[
            Atom("N1", "N", 2, 0, 0),
            Atom("N2", "N", 0, 2, 0),
            Atom("S1", "S", -2, 0, 0),
            Atom("S2", "S", 0, -2, 0),
        ],
    )
    assert site.coordination_number == 4
    assert site.donor_set == "N2S2"


def test_summary():
    site = CoordinationSite(
        metal=Atom("Cu", "Cu", 0, 0, 0),
        ligands=[
            Atom("N1", "N", 2, 0, 0),
            Atom("N2", "N", 0, 2, 0),
        ],
        bonds=[
            Bond("Cu", "N1", 2.0),
            Bond("Cu", "N2", 2.0),
        ],
        title="Test site",
    )
    s = site.summary()
    assert "Test site" in s
    assert "Cu" in s


# ── Geometry classification tests ──────────────────────────────────────────

def _make_tetrahedral():
    """Build an ideal tetrahedral site."""
    ligands = [
        Atom("L1", "N", 1, 1, 1),
        Atom("L2", "N", 1, -1, -1),
        Atom("L3", "N", -1, 1, -1),
        Atom("L4", "N", -1, -1, 1),
    ]
    return CoordinationSite(
        metal=Atom("M", "Cu", 0, 0, 0),
        ligands=ligands,
    )


def _make_square_planar():
    """Build an ideal square planar site."""
    ligands = [
        Atom("L1", "N", 2, 0, 0),
        Atom("L2", "N", 0, 2, 0),
        Atom("L3", "N", -2, 0, 0),
        Atom("L4", "N", 0, -2, 0),
    ]
    return CoordinationSite(
        metal=Atom("M", "Pt", 0, 0, 0),
        ligands=ligands,
    )


def _make_square_pyramidal():
    """Build an ideal square pyramidal site (basal + apical)."""
    ligands = [
        Atom("L1", "N", 2, 0, 0),
        Atom("L2", "N", 0, 2, 0),
        Atom("L3", "N", -2, 0, 0),
        Atom("L4", "N", 0, -2, 0),
        Atom("L5", "N", 0, 0, 2),  # apical
    ]
    return CoordinationSite(
        metal=Atom("M", "Cu", 0, 0, 0),
        ligands=ligands,
    )


def test_tau4_tetrahedral():
    site = _make_tetrahedral()
    t = tau4(site)
    assert t is not None
    assert t > 0.8, f"Expected tau4 > 0.8 for Td, got {t:.3f}"


def test_tau4_square_planar():
    site = _make_square_planar()
    t = tau4(site)
    assert t is not None
    assert t < 0.2, f"Expected tau4 < 0.2 for D4h, got {t:.3f}"


def test_tau5_square_pyramidal():
    site = _make_square_pyramidal()
    t = tau5(site)
    assert t is not None
    assert t < 0.2, f"Expected tau5 ~ 0 for ideal SP, got {t:.3f}"


def test_classify_tetrahedral():
    site = _make_tetrahedral()
    geo = classify(site)
    assert "tetrahedral" in geo.lower()


def test_classify_square_planar():
    site = _make_square_planar()
    geo = classify(site)
    assert "square planar" in geo.lower()

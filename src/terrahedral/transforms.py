"""
Coordinate transforms: fractional↔Cartesian, rotations, projections.
"""

from __future__ import annotations

import math
from typing import Sequence


def frac_to_cart(
    fx: float, fy: float, fz: float,
    a: float, b: float, c: float,
    alpha: float, beta: float, gamma: float,
) -> tuple[float, float, float]:
    """
    Fractional → Cartesian using the standard crystallographic convention.

    Parameters
    ----------
    fx, fy, fz : fractional coordinates
    a, b, c : cell lengths (Å)
    alpha, beta, gamma : cell angles (degrees)

    Returns
    -------
    (x, y, z) in Å
    """
    ar = math.radians(alpha)
    br = math.radians(beta)
    gr = math.radians(gamma)

    cos_a, cos_b, cos_g = math.cos(ar), math.cos(br), math.cos(gr)
    sin_g = math.sin(gr)

    # Guard against degenerate cell (gamma = 0° or 180°)
    if abs(sin_g) < 1e-10:
        raise ValueError(
            f"Degenerate cell angle gamma={gamma}° (sin(gamma)≈0). "
            f"Cannot convert fractional to Cartesian coordinates."
        )

    # Column vectors of the transformation matrix
    ax_c = a
    bx_c = b * cos_g
    by_c = b * sin_g
    cx_c = c * cos_b
    cy_c = c * (cos_a - cos_b * cos_g) / sin_g
    cz_c = c * math.sqrt(
        max(0.0, 1 - cos_a**2 - cos_b**2 - cos_g**2 + 2 * cos_a * cos_b * cos_g)
    ) / sin_g

    x = ax_c * fx + bx_c * fy + cx_c * fz
    y = by_c * fy + cy_c * fz
    z = cz_c * fz

    return (x, y, z)


def distance(p1: Sequence[float], p2: Sequence[float]) -> float:
    """Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))


def angle_at_center(
    p1: Sequence[float],
    center: Sequence[float],
    p2: Sequence[float],
) -> float:
    """Angle (degrees) at *center* between vectors to *p1* and *p2*."""
    v1 = [p1[i] - center[i] for i in range(3)]
    v2 = [p2[i] - center[i] for i in range(3)]
    dot = sum(v1[i] * v2[i] for i in range(3))
    m1 = math.sqrt(sum(x**2 for x in v1))
    m2 = math.sqrt(sum(x**2 for x in v2))
    if m1 * m2 == 0:
        return 0.0
    return math.degrees(math.acos(max(-1.0, min(1.0, dot / (m1 * m2)))))


# ── 3D rotations ────────────────────────────────────────────────────────────

def rotate_y(point: Sequence[float], angle_deg: float) -> tuple[float, float, float]:
    a = math.radians(angle_deg)
    x, y, z = point
    c, s = math.cos(a), math.sin(a)
    return (x * c + z * s, y, -x * s + z * c)


def rotate_x(point: Sequence[float], angle_deg: float) -> tuple[float, float, float]:
    a = math.radians(angle_deg)
    x, y, z = point
    c, s = math.cos(a), math.sin(a)
    return (x, y * c - z * s, y * s + z * c)


def rotate_z(point: Sequence[float], angle_deg: float) -> tuple[float, float, float]:
    a = math.radians(angle_deg)
    x, y, z = point
    c, s = math.cos(a), math.sin(a)
    return (x * c - y * s, x * s + y * c, z)


# ── Projections ─────────────────────────────────────────────────────────────

def project_ortho(
    point_3d: Sequence[float],
    cx: float, cy: float,
    scale: float,
) -> dict:
    """Orthographic projection → {sx, sy, depth}."""
    return {
        "sx": cx + point_3d[0] * scale,
        "sy": cy - point_3d[1] * scale,
        "depth": point_3d[2],
    }


def project_perspective(
    point_3d: Sequence[float],
    cx: float, cy: float,
    scale: float,
    focal: float = 8.0,
) -> dict:
    """Perspective projection → {sx, sy, depth, scale_factor}."""
    d = focal
    s = d / (d - point_3d[2])
    return {
        "sx": cx + point_3d[0] * scale * s,
        "sy": cy - point_3d[1] * scale * s,
        "depth": point_3d[2],
        "scale": s,
    }


def optimal_viewing_angle(
    points: list[Sequence[float]],
    y_range: tuple[int, int] = (0, 360),
    x_range: tuple[int, int] = (-60, 60),
    step: int = 15,
) -> tuple[int, int]:
    """
    Find the Y, X rotation angles that maximize minimum 2D separation
    between projected ligand points.

    Returns (best_y_deg, best_x_deg).
    """
    best_sep = -1.0
    best_angles = (40, -20)

    for yd in range(y_range[0], y_range[1], step):
        for xd in range(x_range[0], x_range[1], step):
            projected = []
            for p in points:
                r = rotate_y(p, yd)
                r = rotate_x(r, xd)
                projected.append((r[0], r[1]))

            min_sep = float("inf")
            for i in range(len(projected)):
                for j in range(i + 1, len(projected)):
                    dx = projected[i][0] - projected[j][0]
                    dy = projected[i][1] - projected[j][1]
                    sep = math.sqrt(dx * dx + dy * dy)
                    min_sep = min(min_sep, sep)

            if min_sep > best_sep:
                best_sep = min_sep
                best_angles = (yd, xd)

    return best_angles

"""
Coordination-shell comparison via Kabsch alignment.

Compare two metal-centred first coordination spheres as 3D point clouds.
After centering on the metal, optimally match donor atoms (Hungarian
assignment), align with Kabsch rotation, and compute RMSD.  Optionally
allow reflection to test pseudo-mirror similarity.

Pure-Python — no numpy/scipy required.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Optional


# ── Linear algebra helpers (3×3, pure Python) ──────────────────────────────

def _mat_mul(A, B):
    """Multiply two 3×3 matrices."""
    return [
        [sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)]
        for i in range(3)
    ]


def _mat_T(A):
    """Transpose a 3×3 matrix."""
    return [[A[j][i] for j in range(3)] for i in range(3)]


def _mat_vec(A, v):
    """Multiply 3×3 matrix by 3-vector."""
    return [sum(A[i][j] * v[j] for j in range(3)) for i in range(3)]


def _det3(M):
    """Determinant of 3×3 matrix."""
    return (
        M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1])
      - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
      + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0])
    )


def _identity():
    return [[1, 0, 0], [0, 1, 0], [0, 0, 1]]


def _jacobi_eigen(M, max_iter=50):
    """
    Eigenvalues and eigenvectors of a 3×3 real symmetric matrix.
    Returns (eigenvalues[3], eigenvectors[3×3]) where eigenvectors
    are columns: V[:,i] is eigenvector for eigenvalue[i].
    Sorted by eigenvalue descending.
    """
    A = [row[:] for row in M]
    V = _identity()

    for _ in range(max_iter):
        # Find largest off-diagonal
        p, q, mx = 0, 1, abs(A[0][1])
        for i, j in [(0, 2), (1, 2)]:
            if abs(A[i][j]) > mx:
                p, q, mx = i, j, abs(A[i][j])
        if mx < 1e-15:
            break

        if abs(A[p][p] - A[q][q]) < 1e-15:
            theta = math.pi / 4
        else:
            theta = 0.5 * math.atan2(2 * A[p][q], A[p][p] - A[q][q])
        c, s = math.cos(theta), math.sin(theta)

        new_A = [row[:] for row in A]
        new_A[p][p] = c*c*A[p][p] + 2*s*c*A[p][q] + s*s*A[q][q]
        new_A[q][q] = s*s*A[p][p] - 2*s*c*A[p][q] + c*c*A[q][q]
        new_A[p][q] = new_A[q][p] = 0.0

        for r in range(3):
            if r == p or r == q:
                continue
            arp = c*A[r][p] + s*A[r][q]
            arq = -s*A[r][p] + c*A[r][q]
            new_A[r][p] = new_A[p][r] = arp
            new_A[r][q] = new_A[q][r] = arq
        A = new_A

        new_V = [row[:] for row in V]
        for r in range(3):
            new_V[r][p] = c*V[r][p] + s*V[r][q]
            new_V[r][q] = -s*V[r][p] + c*V[r][q]
        V = new_V

    eigenvalues = [A[i][i] for i in range(3)]
    order = sorted(range(3), key=lambda i: -eigenvalues[i])
    eigenvalues = [eigenvalues[i] for i in order]
    V_reordered = [[V[r][order[c]] for c in range(3)] for r in range(3)]

    return eigenvalues, V_reordered


def _svd3(H):
    """
    SVD of a 3×3 matrix H = U @ diag(S) @ V^T.
    Returns U (3×3), S (3,), V (3×3).
    Singular values are non-negative.  det(U) and det(V) may be ±1.
    """
    HtH = _mat_mul(_mat_T(H), H)
    eigenvalues, V = _jacobi_eigen(HtH)

    # Singular values = sqrt of eigenvalues of H^T H
    S = [math.sqrt(max(0.0, ev)) for ev in eigenvalues]

    # U = H @ V @ diag(1/S)
    HV = _mat_mul(H, V)
    U = [[0.0] * 3 for _ in range(3)]
    for j in range(3):
        if S[j] > 1e-12:
            for i in range(3):
                U[i][j] = HV[i][j] / S[j]
        else:
            # Zero singular value — set arbitrary orthogonal column
            U[0][j] = 1.0 if j == 0 else 0.0
            U[1][j] = 1.0 if j == 1 else 0.0
            U[2][j] = 1.0 if j == 2 else 0.0

    # Don't force det(U) = +1 here; let Kabsch handle sign correction.
    # But do ensure columns of U are orthonormal via Gram-Schmidt for
    # numerical stability
    for j in range(3):
        # Subtract projection onto previous columns
        for k in range(j):
            dot = sum(U[i][j] * U[i][k] for i in range(3))
            for i in range(3):
                U[i][j] -= dot * U[i][k]
        # Normalize
        norm = math.sqrt(sum(U[i][j]**2 for i in range(3)))
        if norm > 1e-12:
            for i in range(3):
                U[i][j] /= norm

    return U, S, V


# ── Hungarian algorithm ────────────────────────────────────────────────────

def _hungarian(cost_matrix: list[list[float]]) -> list[tuple[int, int]]:
    """
    Solve the assignment problem: find minimum-cost matching.
    cost_matrix is n×m (n ≤ m). Returns list of (row, col) pairs.
    Uses the Jonker-Volgenant / Hungarian algorithm for small matrices.
    """
    n = len(cost_matrix)
    m = len(cost_matrix[0]) if n > 0 else 0
    assert n <= m, "More rows than columns"

    # Pad to square if needed
    size = max(n, m)
    C = [[0.0] * size for _ in range(size)]
    for i in range(n):
        for j in range(m):
            C[i][j] = cost_matrix[i][j]

    # Kuhn-Munkres
    u = [0.0] * (size + 1)
    v = [0.0] * (size + 1)
    p = [0] * (size + 1)
    way = [0] * (size + 1)

    for i in range(1, size + 1):
        p[0] = i
        j0 = 0
        minv = [float('inf')] * (size + 1)
        used = [False] * (size + 1)

        while True:
            used[j0] = True
            i0 = p[j0]
            delta = float('inf')
            j1 = -1

            for j in range(1, size + 1):
                if used[j]:
                    continue
                cur = C[i0 - 1][j - 1] - u[i0] - v[j]
                if cur < minv[j]:
                    minv[j] = cur
                    way[j] = j0
                if minv[j] < delta:
                    delta = minv[j]
                    j1 = j

            for j in range(size + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta

            j0 = j1
            if p[j0] == 0:
                break

        while j0:
            p[j0] = p[way[j0]]
            j0 = way[j0]

    result = []
    for j in range(1, size + 1):
        if p[j] != 0 and p[j] - 1 < n and j - 1 < m:
            result.append((p[j] - 1, j - 1))

    return result


# ── Kabsch alignment ───────────────────────────────────────────────────────

def _kabsch(P: list[list[float]], Q: list[list[float]], allow_reflection: bool = False):
    """
    Compute the optimal rotation (and optionally reflection) that minimises
    RMSD between point sets P and Q (both n×3, already centred).

    Returns (R, rmsd) where R is the 3×3 rotation/reflection matrix.
    """
    n = len(P)
    if n == 0:
        return _identity(), 0.0

    # Cross-covariance matrix H = P^T @ Q
    H = [[0.0] * 3 for _ in range(3)]
    for k in range(n):
        for i in range(3):
            for j in range(3):
                H[i][j] += P[k][i] * Q[k][j]

    U, S, V = _svd3(H)

    # d = sign of det(V @ U^T)
    d = _det3(_mat_mul(V, _mat_T(U)))

    if allow_reflection:
        # Allow improper rotation: R = V @ U^T directly
        # This gives the best possible RMSD including mirror operations
        R = _mat_mul(V, _mat_T(U))
    else:
        # Proper rotation only: correct sign to ensure det(R) = +1
        sign = [1.0, 1.0, 1.0]
        if d < 0:
            sign[2] = -1.0
        D = [[sign[0], 0, 0], [0, sign[1], 0], [0, 0, sign[2]]]
        R = _mat_mul(_mat_mul(V, D), _mat_T(U))

    # Compute RMSD after transformation
    rmsd_sum = 0.0
    for k in range(n):
        rp = _mat_vec(R, P[k])
        for i in range(3):
            rmsd_sum += (rp[i] - Q[k][i]) ** 2
    rmsd = math.sqrt(rmsd_sum / n)

    return R, rmsd


# ── Shell comparison ───────────────────────────────────────────────────────

@dataclass
class ShellComparison:
    """Result of comparing two coordination shells."""
    cn_a: int
    cn_b: int
    donors_a: list[str]
    donors_b: list[str]
    bonds_a: list[float]
    bonds_b: list[float]
    labels_a: list[str]
    labels_b: list[str]
    matched_pairs: list[tuple[int, int]]  # indices into donors_a, donors_b
    rmsd_rotation: float
    rmsd_reflection: float
    percent_improvement: float
    rotation_matrix: list[list[float]]
    reflection_matrix: list[list[float]]
    # Aligned coordinates for visualisation
    coords_a_centred: list[list[float]]   # site A ligands (centred, unrotated)
    coords_b_centred: list[list[float]]   # site B ligands (centred, target)
    coords_a_rotated: list[list[float]]   # site A ligands after Kabsch rotation
    coords_a_reflected: list[list[float]] # site A ligands after Kabsch + reflection
    metal_a: str = ""
    metal_b: str = ""
    element_filter: str = ""              # "" = all, "O" = O-only, etc.


def compare_shells(
    site_a,
    site_b,
    element_filter: str = "",
) -> ShellComparison:
    """
    Compare two CoordinationSite objects.

    Rotation-only: Hungarian match → Kabsch rotation → RMSD.
    Reflection:    Try each axis-flip on site A, re-match and
                   re-align for each, keep the best RMSD.

    Parameters
    ----------
    site_a, site_b : CoordinationSite
    element_filter : if non-empty, only include donors of this element (e.g. "O")
    """
    def _extract(site):
        bd = site.bond_dict
        coords, elements, labels, dists = [], [], [], []
        for lig in site.ligands:
            if element_filter and lig.element != element_filter:
                continue
            coords.append([lig.x - site.metal.x, lig.y - site.metal.y, lig.z - site.metal.z])
            elements.append(lig.element)
            labels.append(lig.residue or lig.label)
            dists.append(bd.get(lig.label, lig.distance_to(site.metal)))
        return coords, elements, labels, dists

    coords_a, elems_a, labels_a, dists_a = _extract(site_a)
    coords_b, elems_b, labels_b, dists_b = _extract(site_b)
    n_a, n_b = len(coords_a), len(coords_b)

    empty = ShellComparison(
        cn_a=n_a, cn_b=n_b, donors_a=elems_a, donors_b=elems_b,
        bonds_a=dists_a, bonds_b=dists_b, labels_a=labels_a, labels_b=labels_b,
        matched_pairs=[], rmsd_rotation=0, rmsd_reflection=0,
        percent_improvement=0,
        rotation_matrix=_identity(), reflection_matrix=_identity(),
        coords_a_centred=coords_a, coords_b_centred=coords_b,
        coords_a_rotated=coords_a, coords_a_reflected=coords_a,
        metal_a=site_a.metal.element, metal_b=site_b.metal.element,
        element_filter=element_filter,
    )
    if n_a == 0 or n_b == 0:
        return empty

    # Ensure P (rows) ≤ Q (cols) for Hungarian
    if n_a <= n_b:
        P_all, Q_all = coords_a, coords_b
        swapped = False
    else:
        P_all, Q_all = coords_b, coords_a
        swapped = True

    def _match_and_align(P_in, Q_in, allow_reflection=False):
        """Run Hungarian + Kabsch on the given point sets."""
        cost = [[sum((P_in[i][k] - Q_in[j][k])**2 for k in range(3))
                 for j in range(len(Q_in))]
                for i in range(len(P_in))]
        pairs = _hungarian(cost)
        P_m = [P_in[r] for r, _ in sorted(pairs)]
        Q_m = [Q_in[c] for _, c in sorted(pairs)]
        R, rmsd = _kabsch(P_m, Q_m, allow_reflection=allow_reflection)
        return R, rmsd, sorted(pairs), P_m, Q_m

    # ── Rotation-only alignment ──
    R_rot, rmsd_rot, pairs_rot, P_matched, Q_matched = _match_and_align(
        P_all, Q_all, allow_reflection=False
    )

    # ── Reflection alignment: try axis-flips on P BEFORE matching ──
    # Try identity + 3 single-axis flips + 3 double-axis flips + full inversion
    # = 8 sign combinations.  Each one reflects P, re-runs Hungarian + Kabsch.
    FLIPS = [
        [1, 1, 1],    # identity (should match rotation result)
        [-1, 1, 1],   # flip X
        [1, -1, 1],   # flip Y
        [1, 1, -1],   # flip Z
        [-1, -1, 1],  # flip XY
        [-1, 1, -1],  # flip XZ
        [1, -1, -1],  # flip YZ
        [-1, -1, -1], # full inversion
    ]

    best_ref_rmsd = float('inf')
    best_ref_R = _identity()
    best_ref_flip = [1, 1, 1]
    best_ref_pairs = pairs_rot

    for flip in FLIPS:
        P_flipped = [[P_all[i][j] * flip[j] for j in range(3)]
                      for i in range(len(P_all))]
        R, rmsd, pairs, _, _ = _match_and_align(P_flipped, Q_all, allow_reflection=False)
        if rmsd < best_ref_rmsd:
            best_ref_rmsd = rmsd
            best_ref_R = R
            best_ref_flip = flip
            best_ref_pairs = pairs

    rmsd_ref = round(best_ref_rmsd, 4)
    rmsd_rot = round(rmsd_rot, 4)

    if rmsd_rot > 1e-9:
        pct = (rmsd_rot - rmsd_ref) / rmsd_rot * 100
    else:
        pct = 0.0

    # Build matched pair indices (map back to original A/B indexing)
    if swapped:
        matched_indices = [(c, r) for r, c in pairs_rot]
    else:
        matched_indices = [(r, c) for r, c in pairs_rot]

    # Compute aligned coordinates for visualisation
    def _apply_all(R, coords):
        return [_mat_vec(R, p) for p in coords]

    src = coords_a if not swapped else coords_b
    coords_rotated = _apply_all(R_rot, src)

    # For reflection: first flip, then rotate
    src_flipped = [[src[i][j] * best_ref_flip[j] for j in range(3)]
                    for i in range(len(src))]
    coords_reflected = _apply_all(best_ref_R, src_flipped)

    return ShellComparison(
        cn_a=n_a, cn_b=n_b,
        donors_a=elems_a, donors_b=elems_b,
        bonds_a=[round(d, 3) for d in dists_a],
        bonds_b=[round(d, 3) for d in dists_b],
        labels_a=labels_a, labels_b=labels_b,
        matched_pairs=matched_indices,
        rmsd_rotation=rmsd_rot,
        rmsd_reflection=rmsd_ref,
        percent_improvement=round(pct, 1),
        rotation_matrix=R_rot,
        reflection_matrix=best_ref_R,
        coords_a_centred=coords_a,
        coords_b_centred=coords_b,
        coords_a_rotated=coords_rotated,
        coords_a_reflected=coords_reflected,
        metal_a=site_a.metal.element,
        metal_b=site_b.metal.element,
        element_filter=element_filter,
    )

"""
Generic Christoffel / Riemann / Ricci computation from a metric.

Used to independently (covariantly) cross-check mini-superspace
Euler-Lagrange results, per the "two independent derivations must
agree" principle in the project plan.
"""
import sympy as sp


def christoffel_symbols(g, ginv, coords):
    n = len(coords)
    Gamma = [[[0] * n for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                s = 0
                for d in range(n):
                    s += ginv[a, d] * (
                        sp.diff(g[d, b], coords[c])
                        + sp.diff(g[d, c], coords[b])
                        - sp.diff(g[b, c], coords[d])
                    )
                Gamma[a][b][c] = sp.simplify(s / 2)
    return Gamma


def riemann_tensor(Gamma, coords):
    """Returns R^a_{b c d} (mixed, first index up)."""
    n = len(coords)
    R = [[[[0] * n for _ in range(n)] for _ in range(n)] for _ in range(n)]
    for a in range(n):
        for b in range(n):
            for c in range(n):
                for d in range(n):
                    term = sp.diff(Gamma[a][b][d], coords[c]) - sp.diff(Gamma[a][b][c], coords[d])
                    for e in range(n):
                        term += Gamma[a][c][e] * Gamma[e][b][d] - Gamma[a][d][e] * Gamma[e][b][c]
                    R[a][b][c][d] = sp.simplify(term)
    return R


def ricci_tensor(Riemann, n):
    """R_{bd} = R^a_{b a d}."""
    Ric = sp.zeros(n, n)
    for b in range(n):
        for d in range(n):
            s = 0
            for a in range(n):
                s += Riemann[a][b][a][d]
            Ric[b, d] = sp.simplify(s)
    return Ric


def ricci_scalar(Ric, ginv, n):
    s = 0
    for a in range(n):
        for b in range(n):
            s += ginv[a, b] * Ric[a, b]
    return sp.simplify(s)


def compute_curvature(g, coords):
    """Convenience wrapper: metric -> (Gamma, Riemann, Ricci_tensor, Ricci_scalar)."""
    n = len(coords)
    ginv = g.inv()
    Gamma = christoffel_symbols(g, ginv, coords)
    Riemann = riemann_tensor(Gamma, coords)
    Ric = ricci_tensor(Riemann, n)
    Rs = ricci_scalar(Ric, ginv, n)
    return Gamma, Riemann, Ric, Rs

"""
H2 (docs/SPEC.md): SECOND, INDEPENDENT derivation route for the GP
background field equations, to cross-check
sympy_layer/proca_minisuperspace.py's route (build a "reduced"
mini-superspace Lagrangian via various integration-by-parts identities,
then take Euler-Lagrange derivatives).

This route builds L_total DIRECTLY from each L_i's own covariant
definition -- NO integration-by-parts shortcut is used anywhere:
  - L2 = G2(X): already non-reduced in both routes (no IBP possible/needed).
  - L3 = G3(X) nabla_mu A^mu: uses div_A literally (proca_minisuperspace.py's
    L3_reduced instead IBPs this into `a^3 A0 G3,X Xdot/N`).
  - L4 = G4(X) R + G4,X[...]: the R-piece uses the Ricci scalar LITERALLY
    (proca_minisuperspace.py's L4_R_reduced instead uses the f(R)-style
    "-6a adot^2 F/N - 6a^2 adot Fdot/N" IBP'd identity); the G4,X-piece
    was already non-reduced/algebraic in both routes.
  - L5 = G5(X) G_{mu nu} nabla^mu A^nu - (1/6)G5,X[...]: the G-piece uses
    the actual Einstein tensor directly (this is now ALSO what
    proca_minisuperspace.py does, after the bug found in its earlier
    Bianchi-identity IBP shortcut -- see build_L5_term()'s docstring --
    so this route re-derives that same piece independently rather than
    importing it, keeping the two routes genuinely separate code paths);
    the cubic-bracket piece was already non-reduced/algebraic in both.

Since this route never drops a total-time-derivative term at all, its
Euler-Lagrange equations are guaranteed (by the standard variational
calculus argument: two Lagrangians differing by a total time derivative
have identical Euler-Lagrange equations for every generalized coordinate)
to match proca_minisuperspace.py's IF AND ONLY IF that file's IBP
reductions are themselves exact identities. Three of those four
reductions (L3, L4-R) were independently verified as EXACT identities
(zero symbolic residual against this direct route) during Stage 2
debugging; the fourth (the old L5 G-piece Bianchi-identity shortcut) was
NOT exact, which is how the bug documented in proca_minisuperspace.py's
build_L5_term() was found in the first place.

Both this file's Euler-Lagrange machinery and proca_minisuperspace.py's
now share the same two corrections (found via this same H2 cross-check
during Stage 2): EL_N needs the full d/dt(dL/dNdot) - dL/dN form (not
naive dL/dN), and EL_a needs the Ostrogradsky-corrected
d/dt(dL/dadot) - dL/da - d^2/dt^2(dL/daddot) form (not the plain
first-order form) -- both because this ansatz's Lagrangian genuinely
depends on Ndot and addot once the true Einstein tensor is used for L5.
"""
import sympy as sp

from proca_minisuperspace import flrw_geometry, covariant_derivative_A, t
from tensor_utils import compute_curvature


def build_lagrangian_direct(G2_expr, G3_expr, G4_expr, G5_expr, d2_val=0):
    """Direct (non-IBP'd) construction of L_total = N*a^3*(L2+L3+L4+L5)
    density, for concrete G2..G5 expressions in Xsym=sp.Symbol('X', positive=True)
    (same interned-symbol convention as proca_minisuperspace.build_lagrangian,
    so the SAME concrete G_i expression objects can be passed to both routes)."""
    coords, g, ginv, Gamma, a_t, N_t = flrw_geometry()
    A0_t = sp.Function('A0')(t)
    A_lower = [A0_t, 0, 0, 0]
    n = 4

    nabla_A, M = covariant_derivative_A(Gamma, ginv, A_lower, coords)
    div_A = sp.simplify(sum(M[i, i] for i in range(n)))
    MM_trace = sp.simplify(sum(M[i, j] * M[j, i] for i in range(n) for j in range(n)))

    def raise1_first(T):
        return sp.Matrix(n, n, lambda a_, b_: sp.simplify(sum(ginv[a_, c] * T[c, b_] for c in range(n))))

    def raise_both(T):
        return sp.Matrix(n, n, lambda a_, b_: sp.simplify(
            sum(ginv[a_, c] * ginv[b_, d] * T[c, d] for c in range(n) for d in range(n))))

    N_up2 = raise_both(nabla_A)
    N_up1 = raise1_first(nabla_A)
    P1 = sp.simplify(sum(nabla_A[r, s] * N_up2[r, s] for r in range(n) for s in range(n)))
    Q1 = sp.simplify(sum(nabla_A[r, s] * N_up2[g_, r] * N_up1[s, g_]
                          for r in range(n) for s in range(n) for g_ in range(n)))
    Q2 = sp.simplify(sum(nabla_A[r, s] * N_up2[g_, r] * M[g_, s]
                          for r in range(n) for s in range(n) for g_ in range(n)))

    Xsym = sp.Symbol('X', positive=True)
    X_t = sp.simplify(-sp.Rational(1, 2) * ginv[0, 0] * A0_t**2)

    def at_bg(expr):
        return expr.subs(Xsym, X_t)

    sqrtmg = N_t * a_t**3

    # --- L2: literal, no IBP possible ---
    L2_density = sqrtmg * at_bg(G2_expr)

    # --- L3: literal div_A, no IBP ---
    L3_density = sqrtmg * at_bg(G3_expr) * div_A

    # --- L4: literal Ricci scalar for the R-piece; algebraic G4,X-piece ---
    _, _, Ric, Rs = compute_curvature(g, coords)
    L4_R_density = sqrtmg * at_bg(G4_expr) * Rs
    L4_deriv_density = sqrtmg * at_bg(sp.diff(G4_expr, Xsym)) * (div_A**2 - MM_trace)

    # --- L5: literal Einstein tensor for the G-piece; algebraic cubic bracket ---
    Einstein_lower = sp.Matrix(n, n, lambda i, j: sp.simplify(Ric[i, j] - sp.Rational(1, 2) * g[i, j] * Rs))
    L5_Gpiece_scalar = sp.simplify(sum(Einstein_lower[mu, nu] * N_up2[mu, nu] for mu in range(n) for nu in range(n)))
    L5_Gpiece_density = sqrtmg * at_bg(G5_expr) * L5_Gpiece_scalar

    d2 = d2_val
    cubic_bracket = (div_A**3 - 3 * d2 * div_A * P1 - 3 * (1 - d2) * div_A * MM_trace
                      + (2 - 3 * d2) * Q1 + 3 * d2 * Q2)
    L5_cubic_density = sqrtmg * (-sp.Rational(1, 6)) * at_bg(sp.diff(G5_expr, Xsym)) * cubic_bracket

    L_total = sp.expand(L2_density + L3_density + L4_R_density + L4_deriv_density
                         + L5_Gpiece_density + L5_cubic_density)

    return {'L_total': L_total, 'a_t': a_t, 'N_t': N_t, 'A0_t': A0_t, 'Xsym': Xsym}


def euler_lagrange_all_direct(data):
    """Same Ostrogradsky/full-EL corrections as
    proca_minisuperspace.euler_lagrange_all -- re-implemented independently
    here (not imported) so this really is a separate code path, not a
    thin wrapper that would silently share a bug with the other route."""
    L = data['L_total']
    a_t, N_t, A0_t = data['a_t'], data['N_t'], data['A0_t']

    Ndot = sp.diff(N_t, t)
    dL_dNdot = sp.diff(L, Ndot)
    EL_N = sp.diff(dL_dNdot, t) - sp.diff(L, N_t)

    adot = sp.diff(a_t, t)
    addot = sp.diff(a_t, t, 2)
    dL_dadot = sp.diff(L, adot)
    dL_da = sp.diff(L, a_t)
    dL_daddot = sp.diff(L, addot)
    EL_a = sp.diff(dL_dadot, t) - dL_da - sp.diff(dL_daddot, t, 2)

    dL_dA0dot = sp.diff(L, sp.diff(A0_t, t))
    dL_dA0 = sp.diff(L, A0_t)
    EL_A0 = sp.diff(dL_dA0dot, t) - dL_dA0

    gauge = {N_t: 1, Ndot: 0, sp.diff(N_t, t, 2): 0}
    EL_N = sp.simplify(EL_N.subs(gauge))
    EL_a = sp.simplify(EL_a.subs(gauge))
    EL_A0 = sp.simplify(EL_A0.subs(gauge))
    return EL_N, EL_a, EL_A0


def run_checks():
    """H2: cross-check this DIRECT route against
    proca_minisuperspace.py's REDUCED (IBP'd) route for several
    independent random G2..G5 choices -- must agree EXACTLY (symbolic
    zero, not just a numeric ratio, since both routes use the identical
    N*a^3-weighted convention with no unknown normalisation between them,
    unlike the H1 comparison against the paper's own be1/be2/be3)."""
    import random
    from proca_minisuperspace import build_lagrangian as build_reduced
    from proca_minisuperspace import euler_lagrange_all as euler_lagrange_all_reduced

    results = []
    Xsym = sp.Symbol('X', positive=True)

    def random_power_law(seed):
        rnd = random.Random(seed)
        b = sp.Rational(rnd.randint(1, 5), rnd.randint(1, 3))
        if rnd.random() < 0.5:
            b = -b
        p = sp.Rational(rnd.randint(1, 3), rnd.randint(1, 2))
        return b * Xsym**p

    for seed in [11, 22]:
        G2c = random_power_law(seed * 7 + 1)
        G3c = random_power_law(seed * 7 + 2)
        G4c = random_power_law(seed * 7 + 3)
        G5c = random_power_law(seed * 7 + 4)

        data_direct = build_lagrangian_direct(G2c, G3c, G4c, G5c, d2_val=0)
        EL_N_d, EL_a_d, EL_A0_d = euler_lagrange_all_direct(data_direct)

        data_reduced = build_reduced(include_L5=True, d2_val=0,
                                      G2_expr=G2c, G3_expr=G3c, G4_expr=G4c, G5_expr=G5c)
        EL_N_r, EL_a_r, EL_A0_r = euler_lagrange_all_reduced(data_reduced)

        ok_N = sp.simplify(EL_N_d - EL_N_r) == 0
        ok_a = sp.simplify(EL_a_d - EL_a_r) == 0
        ok_A0 = sp.simplify(EL_A0_d - EL_A0_r) == 0

        results.append((f'H2 [seed={seed}] EL_N: direct route == reduced route (exact symbolic match)',
                         ok_N, None))
        results.append((f'H2 [seed={seed}] EL_a: direct route == reduced route (exact symbolic match)',
                         ok_a, None))
        results.append((f'H2 [seed={seed}] EL_A0: direct route == reduced route (exact symbolic match)',
                         ok_A0, None))

    return results


if __name__ == '__main__':
    results = run_checks()
    n_pass = 0
    for name, ok, detail in results:
        status = 'PASS' if ok else 'FAIL'
        if ok:
            n_pass += 1
        print(f'[{status}] {name}')
        if not ok:
            print(f'       detail: {detail}')
    print(f'\n{n_pass}/{len(results)} checks passed')

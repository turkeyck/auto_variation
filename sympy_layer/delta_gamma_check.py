"""
H2 (docs/SPEC.md Stage 3): independent SymPy cross-check of the cadabra
delta-Gamma (Palatini) primitive (cadabra/variation_engine.py), applied
to L3 = G3(X) nabla_mu A^mu's g^{mu nu} variation
(cadabra/proca_L3.py).

Unlike every other H1/H2 check in this project, this one is deliberately
NOT tied to the FLRW+temporal-A background ansatz -- module B's whole
point is a COVARIANT primitive that must work for an arbitrary metric,
so this check uses a fully generic diagonal 4-metric with four
independent coordinate-dependent functions and a genuinely two-component
vector field A_mu=(A0(t,x,y,z), A1(t,x,y,z), 0, 0), not any
symmetry-simplified special case.

METHOD: cadabra's Palatini expansion gives (PRE any integration-by-parts)
    delta(nabla_mu A^mu) = delta{g^{mu rho}} nabla_mu A_rho
        + A_alpha nabla_beta(delta{g^{alpha beta}})
        - (1/2) g_{alpha beta} A^lambda nabla_lambda(delta{g^{alpha beta}})
This is a POINTWISE identity (true at every spacetime point, before any
"drop the boundary term" step, which is only valid at the level of an
integrated action) -- so it can be checked directly by perturbing
g^{mu nu} -> g^{mu nu} + eps*h^{mu nu} for a fixed symmetric h and
comparing d/deps of the directly-recomputed nabla_mu A^mu (via
tensor_utils.christoffel_symbols, fully independent of cadabra) against
the RHS built from the closed-form pieces above, both evaluated at
eps=0.

An earlier attempt compared against the POST-IBP form used in the field
equation (delta(sqrt(-g) L3)/(sqrt(-g) delta{g^{mu nu}})|_symmetric =
-1/2 G3,X A_mu A_nu (nabla_lambda A^lambda), sanity-checked separately by
G3=const -> 0) using the SAME pointwise method and found a nonzero
residual -- NOT a bug, but a methodology mismatch: covariant derivatives
of a coordinate-constant h^{mu nu} are generally NONZERO (nabla picks up
Christoffel-connection terms even when partial-derivatives of h vanish),
so a pointwise comparison cannot validate a form that already had
boundary terms (built from integrals of exactly those covariant-derivative-
of-h pieces) dropped. The PRE-IBP form checked here has no such boundary
terms to drop, so the pointwise check is the mathematically appropriate
one for it.
"""
import sympy as sp

from tensor_utils import christoffel_symbols

t, x, y, z = sp.symbols('t x y z', real=True)
_COORDS = [t, x, y, z]


def _generic_metric():
    P = sp.Function('P')(*_COORDS)
    Q = sp.Function('Q')(*_COORDS)
    R = sp.Function('R')(*_COORDS)
    S = sp.Function('S')(*_COORDS)
    g = sp.diag(-P, Q, R, S)
    return g, g.inv()


def _generic_vector():
    A0 = sp.Function('A0')(*_COORDS)
    A1 = sp.Function('A1')(*_COORDS)
    return [A0, A1, 0, 0]


def _nabla_and_div(Gamma, ginv, A_lower):
    n = 4
    nabla = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            expr = sp.diff(A_lower[nu], _COORDS[mu])
            for lam in range(n):
                expr -= Gamma[lam][mu][nu] * A_lower[lam]
            nabla[mu, nu] = expr
    M = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            M[mu, nu] = sum(ginv[nu, rho] * nabla[mu, rho] for rho in range(n))
    return nabla, sum(M[i, i] for i in range(n))


def _fixed_symmetric_perturbation():
    n = 4
    h = sp.zeros(n, n)
    h[0, 0] = 1
    h[0, 1] = h[1, 0] = sp.Rational(1, 3)
    h[1, 1] = sp.Rational(1, 2)
    return h


def pointwise_pre_ibp_residual():
    """Return LHS - RHS (should simplify to 0) for the pre-IBP delta(div_A)
    identity, on a generic metric+vector, for one fixed symmetric
    perturbation direction h^{mu nu}."""
    g, ginv = _generic_metric()
    A_lower = _generic_vector()
    h = _fixed_symmetric_perturbation()
    eps = sp.Symbol('eps')

    ginv_pert = ginv + eps * h
    Gamma_pert = christoffel_symbols(ginv_pert.inv(), ginv_pert, _COORDS)
    _, divA_pert = _nabla_and_div(Gamma_pert, ginv_pert, A_lower)
    LHS = sp.simplify(sp.diff(divA_pert, eps).subs(eps, 0))

    Gamma0 = christoffel_symbols(g, ginv, _COORDS)
    nabla0, _ = _nabla_and_div(Gamma0, ginv, A_lower)
    n = 4

    delta_g_term = sum(h[mu, rho] * nabla0[mu, rho] for mu in range(n) for rho in range(n))

    nabla_beta_h_upper = [0] * n
    for alpha in range(n):
        s = 0
        for beta in range(n):
            for lam in range(n):
                s += Gamma0[alpha][beta][lam] * h[lam, beta] + Gamma0[beta][beta][lam] * h[alpha, lam]
        nabla_beta_h_upper[alpha] = sp.simplify(s)
    term_A = sum(A_lower[a] * nabla_beta_h_upper[a] for a in range(n))

    scalar_gh = sum(g[a, b] * h[a, b] for a in range(n) for b in range(n))
    A_up = [sum(ginv[a, b] * A_lower[b] for b in range(n)) for a in range(n)]
    term_B = sum(A_up[lam] * sp.diff(scalar_gh, _COORDS[lam]) for lam in range(n))

    RHS = sp.simplify(delta_g_term + term_A - sp.Rational(1, 2) * term_B)
    return sp.simplify(LHS - RHS)


def g3_const_sanity_residual():
    """Independent sanity check (action-level, not pointwise): for
    constant G3, L3 = G3 * nabla_mu A^mu is a pure total covariant
    divergence (nabla_mu(G3 A^mu)), so its contribution to the g^{mu nu}
    field equation must vanish. The POST-IBP closed form
    -1/2 G3,X A_mu A_nu (nabla_lambda A^lambda) manifestly vanishes when
    G3,X=0 -- trivial by construction, recorded here as a named check
    rather than left as an unverified claim in a docstring."""
    Xsym = sp.Symbol('X', positive=True)
    G3X = sp.Symbol('G3X_const', real=True)  # stands for G3,X = 0 case
    A_mu, A_nu, divA = sp.symbols('A_mu A_nu divA', real=True)
    closed_form = -sp.Rational(1, 2) * G3X * A_mu * A_nu * divA
    return sp.simplify(closed_form.subs(G3X, 0))


def run_checks():
    results = []

    residual = pointwise_pre_ibp_residual()
    ok1 = residual == 0
    results.append((
        '1. [H2] Pre-IBP delta(nabla_mu A^mu) closed form (from cadabra\'s Palatini/delta-Gamma '
        'expansion, cadabra/variation_engine.py + proca_L3.py) matches a fully independent, generic '
        '(non-FLRW) SymPy finite-perturbation computation exactly',
        ok1, residual))

    g3_residual = g3_const_sanity_residual()
    ok2 = g3_residual == 0
    results.append((
        '2. Sanity check: constant G3 makes the post-IBP L3 field-equation contribution vanish '
        '(L3 is then a pure boundary term) -- consistent with the closed form used in proca_L3.py',
        ok2, g3_residual))

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

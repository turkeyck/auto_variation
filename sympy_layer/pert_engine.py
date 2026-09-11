"""
Perturbation engine primitives (plan section "引擎必備的四個零件"):

  1. eps_truncate       - drop O(eps^n) and higher from a bookkeeping expansion
  2. FOURIER_RULES       - w0/w1/w2 pairing rules for a mode function W(x)
                           with w0=W, w1=(one spatial derivative), w2=(Laplacian)
  3. background_reduce   - recursively rewrite addot -> a(H^2+Hdot) (and any
                           higher a-derivative reached by repeated d/dt of
                           that identity)
  4. euler_lagrange_1d   - Euler-Lagrange derivative for a 1-field point
                           Lagrangian, and kinetic_matrix() to read off the
                           quadratic-action kinetic coefficient K = a^3 c_s^2
                           normalisation from it.

Demonstrated on Maldacena's single canonical scalar field result:
    L_2 = a^3 eps_H [zetadot^2 - (partial zeta)^2 / a^2]     (Mpl = 1)
  ->  K = a^3 eps_H ,  c_s^2 = 1
"""
import sympy as sp

t = sp.Symbol('t', real=True)
k = sp.Symbol('k', positive=True)


# ---------- 1. epsilon-truncation ----------
def eps_truncate(expr, eps, order):
    """Keep terms up to and including eps**(order-1); drop eps**order and higher."""
    poly = sp.Poly(sp.expand(expr), eps) if expr.has(eps) else None
    if poly is None:
        return expr
    terms = poly.all_terms()
    deg = poly.degree()
    out = 0
    for monom, coeff in terms:
        power = monom[0]
        if power < order:
            out += coeff * eps**power
    return out


# ---------- 2. Fourier pairing rules ----------
w0, w1, w2 = sp.symbols('w0 w1 w2', real=True)
FOURIER_RULES = {
    w0**2: sp.Integer(1),
    w1**2: k**2,
    w0 * w2: -k**2,
    w2 * w0: -k**2,
    w0 * w1: sp.Integer(0),
    w1 * w0: sp.Integer(0),
    w1 * w2: sp.Integer(0),
    w2 * w1: sp.Integer(0),
}


def fourier_reduce(expr):
    return sp.expand(expr).subs(FOURIER_RULES)


# ---------- 3. background-equation reduction ----------
def background_reduce(expr, a, H, Hdot_expr=None, max_order=6):
    """Recursively rewrite d^n a/dt^n (n>=2) using addot = a(H^2+Hdot),
    Hdot = dH/dt (kept symbolic unless Hdot_expr supplied)."""
    a_t = a
    H_sym = H
    addot_identity = a_t * (H_sym**2 + sp.Derivative(H_sym, t))
    e = expr
    # substitute highest derivatives first so repeated differentiation of the
    # identity itself eliminates order-by-order
    for n in range(max_order, 1, -1):
        target = sp.Derivative(a_t, (t, n))
        if e.has(target):
            # d^n a/dt^n = d^{n-2}/dt^{n-2} [a(H^2+Hdot)]
            repl = sp.diff(addot_identity, t, n - 2)
            e = e.subs(target, repl)
    e = sp.expand(e.doit())
    if Hdot_expr is not None:
        e = e.subs(sp.Derivative(H_sym, t), Hdot_expr)
    return e


# ---------- 4. Euler-Lagrange + kinetic matrix ----------
def euler_lagrange_1d(L, q_t, tvar):
    dL_dqdot = sp.diff(L, sp.diff(q_t, tvar))
    dL_dq = sp.diff(L, q_t)
    return sp.diff(dL_dqdot, tvar) - dL_dq


def kinetic_coefficient(L, q_t, tvar):
    """K = (1/2) d^2 L / d(qdot)^2 -- the coefficient of qdot^2 in a
    Lagrangian quadratic in qdot."""
    qdot = sp.diff(q_t, tvar)
    return sp.simplify(sp.diff(L, qdot, 2) / 2)


def gradient_coefficient(L, q_t, tvar, kvar):
    """Coefficient of -k^2 q^2 in a Fourier-reduced mode Lagrangian
    (i.e. -(1/2) dL/d(k^2) at fixed q, assuming L is linear in k^2)."""
    q = q_t
    return sp.simplify(-sp.diff(L, kvar**2))


def run_checks():
    results = []

    # ---- Check 1: background_reduce correctly eliminates addot via
    # the Friedmann/Hubble identity, verified against direct substitution
    a_t = sp.Function('a')(t)
    H_t = sp.Function('H')(t)
    addot = sp.diff(a_t, t, 2)
    expr = 3 * addot / a_t + 1  # e.g. would appear in a Raychaudhuri-like combo
    reduced = background_reduce(expr, a_t, H_t)
    expected = 3 * (H_t**2 + sp.diff(H_t, t)) + 1
    ok1 = sp.simplify(reduced - expected) == 0
    results.append(('1. background_reduce: addot -> a(H^2+Hdot) applied correctly', ok1, reduced))

    # ---- Check 2: Fourier pairing engine reproduces the standard
    # (partial_i zeta)(partial^i zeta) -> k^2 zeta_k^2 rule used throughout
    # cosmological perturbation theory (plan "Fourier化簡" component)
    zeta_k = sp.Symbol('zeta_k', real=True)
    zetadot_k = sp.Symbol('zetadot_k', real=True)
    density_grad2 = (zeta_k * w1)**2               # (partial zeta)^2 pattern
    density_zdot2 = (zetadot_k * w0)**2
    reduced_grad = fourier_reduce(density_grad2)
    reduced_zdot = fourier_reduce(density_zdot2)
    ok2 = (sp.simplify(reduced_grad - k**2 * zeta_k**2) == 0
           and sp.simplify(reduced_zdot - zetadot_k**2) == 0)
    results.append(('2. Fourier rules: (partial zeta)^2 -> k^2 zeta_k^2, zetadot^2 -> zetadot_k^2',
                     ok2, (reduced_grad, reduced_zdot)))

    # ---- Check 3 (plan 3.f, Maldacena): build the mode Lagrangian for the
    # well-known canonical single-field quadratic curvature-perturbation
    # action L_2 = a^3 eps_H [zetadot^2 - (partial zeta)^2/a^2], run it
    # through the EL + kinetic-coefficient machinery, and confirm
    # K = a^3 eps_H, c_s^2 = 1.
    a_sym = sp.Function('a')(t)
    epsH = sp.Symbol('epsilon_H', positive=True)
    zk = sp.Function('zeta_k')(t)
    zkdot = sp.diff(zk, t)

    L_mode = a_sym**3 * epsH * (zkdot**2 - k**2 / a_sym**2 * zk**2)

    K = kinetic_coefficient(L_mode, zk, t)
    # L_mode = a^3 epsH zkdot^2 - epsH a k^2 zk^2  -> coefficient of k^2 zk^2 is -epsH*a
    coeff_k2 = sp.simplify(sp.diff(L_mode, k, 2) / 2)  # d^2L/dk^2 /2 = coefficient of k^2*zeta^2 term
    # canonical form L = K zkdot^2 - (K c_s^2 / a^2) k^2 zk^2  ->  c_s^2 = -coeff_k2 * a^2 / (K * zk^2)
    cs2 = sp.simplify((-coeff_k2 / zk**2) * a_sym**2 / K)

    EOM = euler_lagrange_1d(L_mode, zk, t)
    EOM_expected = sp.diff(2 * a_sym**3 * epsH * zkdot, t) + 2 * epsH * a_sym * k**2 * zk
    ok3a = sp.simplify(K - a_sym**3 * epsH) == 0
    ok3b = sp.simplify(cs2 - 1) == 0
    ok3c = sp.simplify(EOM - EOM_expected) == 0
    ok3 = ok3a and ok3b and ok3c
    results.append(('3. Maldacena: K = a^3 eps_H, c_s^2 = 1, mode EOM matches direct EL', ok3,
                     (K, cs2, sp.simplify(EOM - EOM_expected))))

    # ---- Check 4: eps_truncate drops higher-order bookkeeping terms
    eps = sp.Symbol('eps')
    x = sp.Symbol('x')
    expr4 = x + eps * x**2 + eps**2 * x**3 + eps**3 * x**4
    trunc4 = eps_truncate(expr4, eps, 3)  # keep up to eps^2
    ok4 = sp.simplify(trunc4 - (x + eps * x**2 + eps**2 * x**3)) == 0
    results.append(('4. eps_truncate drops O(eps^3) and higher', ok4, trunc4))

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

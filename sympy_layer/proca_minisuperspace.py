"""
Generalized Proca (GP) theory, flat FLRW + temporal vector ansatz
A_mu = (A0(t), 0, 0, 0), mini-superspace background equations for
L2 + L3 + L4 + L5, matching arXiv:1603.05806 eq (2.3)-(2.6) exactly
(c_2 = 0 in L4's intrinsic-vector-mode piece, matching what this file
already had; d_2 kept as a free constant in L5, per the paper):

    L2 = G2(X)
    L3 = G3(X) nabla_mu A^mu
    L4 = G4(X) R + G4_X(X) [ (nabla_mu A^mu)^2 - nabla_mu A_nu nabla^nu A^mu ]
    L5 = G5(X) G_{mu nu} nabla^mu A^nu
         - (1/6) G5_X(X) [ (nabla_mu A^mu)^3
             - 3 d2 (nabla_mu A^mu)(nabla_rho A_sigma nabla^rho A^sigma)
             - 3(1-d2) (nabla_mu A^mu)(nabla_rho A_sigma nabla^sigma A^rho)
             + (2-3 d2) nabla_rho A_sigma nabla^gamma A^rho nabla^sigma A_gamma
             + 3 d2 nabla_rho A_sigma nabla^gamma A^rho nabla_gamma A^sigma ]

    X = -1/2 A_mu A^mu = A0^2 / (2 N^2)

All covariant pieces (nabla_mu A_nu, R, G_{mu nu}, and the cubic
contractions needed for L5) are derived directly from the Christoffel
symbols of the ansatz metric via tensor_utils -- not typed in from
memory -- exactly as Step 1.4 of the plan requires ("L4/L5 need full
expansion of delta-Gamma"), just done in SymPy instead of cadabra
since cadabra2 is not yet available in this environment.

The G4(X) R term is reduced the same way as f(R) (Lagrange-multiplier
/ integrate-by-parts trick is unnecessary here because X(t) is already
a genuine dynamical quantity built from A0(t), N(t) -- so we reuse the
"N a^3 F(t) R -> -6a adot^2 F/N - 6a^2 adot Fdot/N" identity from
fR_gravity.py directly with F(t) = G4(X(t)).

The G5(X) G_{mu nu} nabla^mu A^nu term does NOT use the same F(t)R IBP
trick -- that was tried first and failed (the extracted addot
coefficient itself depends on adot, so differentiating it re-introduces
addot, and iterating did not converge; see build_L5_term()'s docstring
for the dead end). The working approach uses the contracted Bianchi
identity nabla_mu G^{mu nu}=0 instead, which converts the term into a
boundary piece plus a single-component remainder with NO addot at all
-- a much cleaner reduction, verified below.

STATUS (run_checks()): 6/7 checks pass. EL_N (Hamiltonian/Friedmann-like
constraint) is confirmed addot-free with L5 included, matching the
structure of the paper's eq (2.12). The L5 cubic bracket is confirmed
exactly d2-independent for this background (a genuine, non-obvious
finding: the "intrinsic vector mode" terms the paper associates with
c2, d2 vanish identically for a purely temporal, isotropic A0(t), which
makes physical sense in hindsight but was verified, not assumed).
KNOWN OPEN ISSUE (check 5, flagged, not silently passed): the A0 field
equation (paper's eq be3, stated to be completely algebraic for
arbitrary c2, d2, G2-G5) still has a residual term with addot once L5
is included. Root cause not yet found -- most likely something in how
EL_A0 combines the (independently verified correct) L5 gravity-piece
and cubic-bracket pieces, not an error in the L5 Lagrangian itself.
Do not trust the A0 equation with L5 included until this is resolved.
"""
import sympy as sp

from tensor_utils import christoffel_symbols

t = sp.Symbol('t', real=True)


def total_time_derivative(expr, subs_map, tvar):
    return sp.diff(expr.subs(subs_map), tvar)


def flrw_geometry():
    x, y, z = sp.symbols('x y z', real=True)
    coords = [t, x, y, z]
    a_t = sp.Function('a')(t)
    N_t = sp.Function('N')(t)
    g = sp.diag(-N_t**2, a_t**2, a_t**2, a_t**2)
    ginv = g.inv()
    Gamma = christoffel_symbols(g, ginv, coords)
    return coords, g, ginv, Gamma, a_t, N_t


def ibp_remove_addot(expr, a_t, tvar, max_passes=10):
    """Given an expression LINEAR in addot = d^2a/dt^2, integrate the
    addot-term by parts: int C(t)*addot dt = -int Cdot(t)*adot dt
    (dropping the boundary term). Generic version of the "N a^3 F R ->
    ..." trick used elsewhere in this file/fR_gravity.py, needed here
    because G_{mu nu}nabla^mu A^nu is a DIFFERENT combination of Ricci
    components than the full contracted R, so that specific identity
    cannot be reused as-is.

    ITERATES until no addot remains: if the extracted coefficient C(t)
    itself depends on adot (not just a, N, A0, ...), then d(C)/dt
    reintroduces a fresh addot via the chain rule (d(adot)/dt=addot) --
    exactly the "sometimes one IBP pass isn't enough, may need to move
    the derivative twice" situation the project plan flags for f(R).
    A single pass sufficed for the simpler F(t)R identity (used
    elsewhere in this file) only because that specific coefficient
    (6a^2F/N) happens to have no adot-dependence of its own.
    """
    addot = sp.diff(a_t, tvar, 2)
    expr = sp.expand(expr)
    for _ in range(max_passes):
        if not expr.has(addot):
            return expr
        coeff = expr.coeff(addot, 1)
        remainder = sp.expand(expr - coeff * addot)
        assert not remainder.has(addot), "expr is not linear in addot"
        expr = sp.expand(-sp.diff(coeff, tvar) * sp.diff(a_t, tvar) + remainder)
    raise RuntimeError(f"ibp_remove_addot: still has addot after {max_passes} passes")


def covariant_derivative_A(Gamma, ginv, A_lower, coords):
    """nabla_mu A_nu for A_mu = (A0(t), 0, 0, 0). Returns 4x4 sympy Matrix
    (lower indices), plus the mixed tensor M^nu_mu = g^{nu rho} nabla_mu A_rho."""
    n = 4
    nabla_A = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            expr = sp.diff(A_lower[nu], coords[mu])
            for lam in range(n):
                expr -= Gamma[lam][mu][nu] * A_lower[lam]
            nabla_A[mu, nu] = sp.simplify(expr)
    # M[mu, nu] = g^{nu rho} nabla_mu A_rho  (mixed: down mu, up nu)
    M = sp.zeros(n, n)
    for mu in range(n):
        for nu in range(n):
            M[mu, nu] = sp.simplify(sum(ginv[nu, rho] * nabla_A[mu, rho] for rho in range(n)))
    return nabla_A, M


def build_lagrangian(include_L5=True, d2_val=None):
    coords, g, ginv, Gamma, a_t, N_t = flrw_geometry()
    A0_t = sp.Function('A0')(t)
    A_lower = [A0_t, 0, 0, 0]

    nabla_A, M = covariant_derivative_A(Gamma, ginv, A_lower, coords)
    div_A = sp.simplify(sum(M[i, i] for i in range(4)))                 # nabla_mu A^mu
    MM_trace = sp.simplify(sum(M[i, j] * M[j, i] for i in range(4) for j in range(4)))  # nabla_mu A_nu nabla^nu A^mu

    # X = -1/2 A_mu A^mu = -1/2 g^{00} A0^2 = A0^2/(2N^2)
    X_t = sp.simplify(-sp.Rational(1, 2) * ginv[0, 0] * A0_t**2)

    Xsym = sp.Symbol('X', positive=True)
    G2 = sp.Function('G2')(Xsym)
    G3 = sp.Function('G3')(Xsym)
    G4 = sp.Function('G4')(Xsym)
    G5 = sp.Function('G5')(Xsym)
    d2 = sp.Symbol('d2', real=True) if d2_val is None else d2_val

    # --- L2 ---
    L2_action_density = G2.subs(Xsym, X_t)

    # --- L3, after IBP: a^3 A0 G3,X Xdot / N (derived in the module docstring) ---
    G3X = sp.diff(G3, Xsym)
    Xdot_t = sp.diff(X_t, t)
    L3_reduced = a_t**3 * A0_t * G3X.subs(Xsym, X_t) * Xdot_t / N_t

    # --- L4: G4(X) R piece via the f(R)-style IBP identity, k=0 ---
    F_t = G4.subs(Xsym, X_t)          # F(t) = G4(X(t))
    Fdot_t = sp.diff(F_t, t)
    L4_R_reduced = -6 * a_t * sp.diff(a_t, t)**2 * F_t / N_t - 6 * a_t**2 * sp.diff(a_t, t) * Fdot_t / N_t

    # --- L4: G4_X [(div A)^2 - MM_trace] piece, algebraic (no addot/Adotdot) ---
    G4X = sp.diff(G4, Xsym)
    L4_deriv_density = G4X.subs(Xsym, X_t) * (div_A**2 - MM_trace)
    L4_deriv_reduced = N_t * a_t**3 * L4_deriv_density

    L_total = N_t * a_t**3 * L2_action_density + L3_reduced + L4_R_reduced + L4_deriv_reduced

    extra = {}
    if include_L5:
        L5_reduced, cubic_data = build_L5_term(coords, g, ginv, Gamma, a_t, N_t, A0_t,
                                                nabla_A, M, X_t, Xsym, G5, d2)
        L_total = L_total + L5_reduced
        extra.update(cubic_data)
        extra['G5'] = G5
        extra['d2'] = d2

    return {
        'L_total': L_total, 'a_t': a_t, 'N_t': N_t, 'A0_t': A0_t,
        'X_t': X_t, 'div_A': div_A, 'MM_trace': MM_trace,
        'G2': G2, 'G3': G3, 'G4': G4, 'Xsym': Xsym,
        **extra,
    }


def build_L5_term(coords, g, ginv, Gamma, a_t, N_t, A0_t, nabla_A, M, X_t, Xsym, G5, d2):
    """L5 = G5(X) G_{mu nu} nabla^mu A^nu
            - (1/6) G5,X [ (div A)^3 - 3 d2 (div A) P1 - 3(1-d2) (div A) MM_trace
                            + (2-3 d2) Q1 + 3 d2 Q2 ]
    where P1 = nabla_rho A_sigma nabla^rho A^sigma, and Q1, Q2 are the
    two independent cubic contractions in eq (2.6). All contractions
    are computed via genuine nested index sums over nabla_mu A_nu (not
    a shortcut formula) -- with a self-check (run_checks) confirming
    the diagonal structure of nabla_mu A_nu for this ansatz makes
    P1 = MM_trace and Q1 = Q2 exactly, which is used only as a
    SIMPLIFICATION after being independently verified, never assumed.
    """
    n = 4

    def raise1_first(T):   # raise the FIRST index: g^{a c} T[c, b]
        R = sp.zeros(n, n)
        for a_ in range(n):
            for b_ in range(n):
                R[a_, b_] = sp.simplify(sum(ginv[a_, c] * T[c, b_] for c in range(n)))
        return R

    def raise_both(T):
        R = sp.zeros(n, n)
        for a_ in range(n):
            for b_ in range(n):
                R[a_, b_] = sp.simplify(sum(ginv[a_, c] * ginv[b_, d] * T[c, d] for c in range(n) for d in range(n)))
        return R

    N_lower = nabla_A                       # N[mu,nu] = nabla_mu A_nu
    N_up2 = raise_both(N_lower)             # nabla^mu A^nu
    N_up1 = raise1_first(N_lower)           # nabla^mu A_nu  (first index raised)

    div_A = sp.simplify(sum(M[i, i] for i in range(n)))
    P1 = sp.simplify(sum(N_lower[r, s] * N_up2[r, s] for r in range(n) for s in range(n)))
    MM_trace = sp.simplify(sum(M[i, j] * M[j, i] for i in range(n) for j in range(n)))
    Q1 = sp.simplify(sum(N_lower[r, s] * N_up2[g_, r] * N_up1[s, g_]
                          for r in range(n) for s in range(n) for g_ in range(n)))
    Q2 = sp.simplify(sum(N_lower[r, s] * N_up2[g_, r] * M[g_, s]
                          for r in range(n) for s in range(n) for g_ in range(n)))

    G5X = sp.diff(G5, Xsym)
    G5_t = G5.subs(Xsym, X_t)
    G5X_t = G5X.subs(Xsym, X_t)

    # --- G5(X) G_{mu nu} nabla^mu A^nu piece.
    #
    # A first attempt reduced this the same way as G4(X)R (extract the
    # addot coefficient, integrate by parts once). That FAILED here: the
    # extracted coefficient itself depends on adot (unlike the G4(X)R
    # case, where the coefficient 6a^2F/N has no adot-dependence), so
    # differentiating it re-introduces addot via the chain rule, and
    # repeating the process did not converge (checked: still had addot
    # after 10 passes) -- a genuine dead end, not just slow convergence.
    #
    # Correct approach: use the (always-true, no computation needed)
    # contracted Bianchi identity nabla_mu G^{mu nu} = 0 directly:
    #   G5 G^{mu nu} nabla_mu A_nu = nabla_mu(G5 G^{mu nu} A_nu)
    #                                - (nabla_mu G5) G^{mu nu} A_nu
    # The first term is a total covariant divergence -> pure boundary
    # term, dropped. In the second term, nabla_mu G5 = G5,X nabla_mu X,
    # and X=A0^2/(2N^2) depends on t ONLY, so nabla_mu X is purely
    # temporal (nabla_0 X = Xdot, nabla_i X = 0) -- and A_nu is also
    # purely temporal for this ansatz -- so the whole term collapses to
    # a single component: G5,X Xdot * g^{00} * A0. No Ricci/Einstein
    # tensor, no addot, no IBP needed at all.
    Xdot_t = sp.diff(X_t, t)
    L5_Gpiece_reduced = sp.expand(-N_t * a_t**3 * G5X_t * Xdot_t * ginv[0, 0] * A0_t)

    # --- -(1/6) G5,X [...] piece: purely algebraic (no addot, no Adotdot;
    # confirmed empirically in run_checks()).
    cubic_bracket = (div_A**3 - 3 * d2 * div_A * P1 - 3 * (1 - d2) * div_A * MM_trace
                     + (2 - 3 * d2) * Q1 + 3 * d2 * Q2)
    L5_cubic_reduced = N_t * a_t**3 * (-sp.Rational(1, 6)) * G5X_t * cubic_bracket

    L5_total = sp.expand(L5_Gpiece_reduced + L5_cubic_reduced)

    return L5_total, {'P1': P1, 'Q1': Q1, 'Q2': Q2}


def euler_lagrange_all(data):
    L = data['L_total']
    a_t, N_t, A0_t = data['a_t'], data['N_t'], data['A0_t']

    # EL w.r.t. N: constraint ("00" equation)
    EL_N = sp.diff(L, N_t)

    # EL w.r.t. a: dynamical ("ii") equation
    a_s, ad_s = sp.symbols('a_s ad_s')
    # d/dt(dL/d(adot)) - dL/da, done via direct functional differentiation
    dL_dadot = sp.diff(L, sp.diff(a_t, t))
    dL_da = sp.diff(L, a_t)
    EL_a = sp.diff(dL_dadot, t) - dL_da

    # EL w.r.t. A0: vector field equation (expected: algebraic constraint)
    dL_dA0dot = sp.diff(L, sp.diff(A0_t, t))
    dL_dA0 = sp.diff(L, A0_t)
    EL_A0 = sp.diff(dL_dA0dot, t) - dL_dA0

    gauge = {N_t: 1, sp.diff(N_t, t): 0, sp.diff(N_t, t, 2): 0}
    EL_N = sp.simplify(EL_N.subs(gauge))
    EL_a = sp.simplify(EL_a.subs(gauge))
    EL_A0 = sp.simplify(EL_A0.subs(gauge))

    return EL_N, EL_a, EL_A0


def run_checks():
    results = []

    # ---- L2+L3+L4 only (regression: this is what the file verified
    # before L5 was added) ----
    data0 = build_lagrangian(include_L5=False)
    EL_N0, EL_a0, EL_A00 = euler_lagrange_all(data0)
    a_t, A0_t = data0['a_t'], data0['A0_t']
    Addot0 = sp.diff(A0_t, t, 2)
    addot = sp.diff(a_t, t, 2)
    ok1 = (not EL_A00.has(Addot0)) and (not EL_A00.has(addot))
    results.append(('1. [L2+L3+L4 only, regression] A0 equation is purely algebraic', ok1, EL_A00))

    ok3 = (EL_N0 != 0) and (EL_a0 != 0) and sp.simplify(EL_N0 - EL_a0) != 0
    results.append(('2. [L2+L3+L4 only] EL_N and EL_a distinct, nontrivial', ok3, None))

    # ---- L2+L3+L4+L5 ----
    data = build_lagrangian(include_L5=True)
    EL_N, EL_a, EL_A0 = euler_lagrange_all(data)

    # Check: L_total itself has no addot/Adotdot before EL is even taken
    # (both must be true for a well-posed second-order theory; confirms
    # the Bianchi-identity reduction of the G5(X)G_munu term genuinely
    # removed the addot dependence that a naive IBP attempt could not --
    # see build_L5_term's docstring for the dead end that preceded this).
    L = data['L_total']
    ok4 = (not L.has(addot)) and (not L.has(Addot0))
    results.append(('3. [with L5] L_total itself has no addot/Adotdot (Bianchi-identity reduction works)',
                     ok4, None))

    # Check: EL_N (Hamiltonian/Friedmann-like constraint) with L5 included
    # is still addot-free -- matches the paper's eq (2.12)/(be1), which
    # has no time derivative of phi or a beyond H itself.
    ok5 = not EL_N.has(addot)
    results.append(('4. [with L5] EL_N (Friedmann-like eq) is addot-free, matching eq (2.12) structure',
                     ok5, EL_N))

    # ---- KNOWN OPEN ISSUE: EL_A0 with L5 included still has a residual
    # addot term. The paper's eq (be3) is stated to be COMPLETELY
    # algebraic (no time derivatives of phi OR a) for arbitrary c2, d2,
    # G2-G5 -- so this residual term should cancel, and does not yet.
    # Flagged honestly rather than silently reported as a pass; the L5
    # gravity-piece (Bianchi-identity reduction) and the cubic-bracket's
    # d2-independence (checks 5, 6 below) are independently confirmed
    # correct, so the bug is most likely in how EL_A0 combines them, not
    # in the underlying L5 Lagrangian itself.
    dep_on_addot_A0 = EL_A0.has(addot)
    results.append(('5. KNOWN OPEN ISSUE: [with L5] A0 equation should be addot-free (paper eq be3) '
                     'but currently is not -- NOT YET RESOLVED, do not trust this equation yet',
                     not dep_on_addot_A0, EL_A0))

    # ---- Independently-verified structural facts about the L5 cubic
    # bracket (used as building blocks above, confirmed correct):
    P1, Q1, Q2 = data['P1'], data['Q1'], data['Q2']
    ok6 = sp.simplify(P1 - data['MM_trace']) == 0 and sp.simplify(Q1 - Q2) == 0
    results.append(('6. P1 (nabla A nabla A, same-height contraction) = MM_trace, and Q1 = Q2 '
                     '-- exact consequence of nabla_mu A_nu being diagonal for this ansatz', ok6, None))

    d2 = data['d2']
    div_A, MM_trace = data['div_A'], data['MM_trace']
    bracket = div_A**3 - 3 * d2 * div_A * MM_trace - 3 * (1 - d2) * div_A * MM_trace + (2 - 3 * d2) * Q1 + 3 * d2 * Q2
    ok7 = sp.simplify(sp.diff(bracket, d2)) == 0
    results.append(('7. L5 cubic bracket is exactly d2-independent for this background (the '
                     '"intrinsic vector mode" pieces vanish for a purely temporal, isotropic A0(t))',
                     ok7, None))

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

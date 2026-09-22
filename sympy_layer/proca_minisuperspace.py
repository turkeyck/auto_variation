"""
Generalized Proca (GP) theory, flat FLRW + temporal vector ansatz
A_mu = (A0(t), 0, 0, 0) [LOWER index -- see the phi/A0 convention note
below], mini-superspace background equations for L2 + L3 + L4 + L5. The
L2-L5 functional forms below coincide, on this background ansatz, with
BOTH arXiv:1603.05806 eq (2.3)-(2.6) at c2=0 (L4) and d2 arbitrary but
proven background-irrelevant (L5, see check 7) AND arXiv:1703.09573's
own L2-L5 (eq. L2-L5 in that paper's Sec. II) minus the g5(X) term, which
that paper's own text states does not contribute to the background
(astrophv2.tex lines 323-326, quoted in
verification/paper_equations/arxiv_1703_09573.py):

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
memory -- done in SymPy (component-level, final-ansatz substitution;
see docs/SPEC.md module J) rather than an abstract-index covariant engine
(docs/SPEC.md module B, not yet built).

PHI / A0 INDEX-CONVENTION NOTE: this file's A0(t) is the LOWER-index
component A_0 (required by covariant_derivative_A()'s use of the
standard nabla_mu A_nu = partial_mu A_nu - Gamma^lam_mu_nu A_lam formula,
which needs lower components as input). arXiv:1703.09573 states its own
ansatz as A^mu=(phi(t),0,0,0) -- UPPER index. At N=1, phi_paper = A^0 =
g^{00} A_0 = -A_0. Any code comparing this file's background equations
against that paper's be1/be2/be3 must pass -A0_t as phi, not A0_t
directly (verification/paper_equations/arxiv_1703_09573.py's docstring
and sympy_layer/proca_background_1703_09573.py do this correctly; an
earlier draft of that comparison script did not, and appeared to show a
completely different functional dependence for the G3 sector until this
was found and fixed -- even powers of A0/phi, e.g. everything built from
X, are unaffected by the sign, which is why the bug went unnoticed until
an isolated single-term test exposed it).

The G4(X) R term is reduced the same way as f(R) (Lagrange-multiplier
/ integrate-by-parts trick is unnecessary here because X(t) is already
a genuine dynamical quantity built from A0(t), N(t) -- so we reuse the
"N a^3 F(t) R -> -6a adot^2 F/N - 6a^2 adot Fdot/N" identity from
fR_gravity.py directly with F(t) = G4(X(t)).

The G5(X) G_{mu nu} nabla^mu A^nu term: TWO REDUCTION ATTEMPTS FAILED
before the current, verified approach (see build_L5_term()'s docstring
for full detail). (1) The F(t)R-style IBP trick failed outright (does
not converge). (2) A contracted-Bianchi-identity IBP trick appeared to
work (passed this file's own internal regression for a long time) but
was PROVEN WRONG by an H1 comparison against arXiv:1703.09573's own be1
(docs/SPEC.md Stage 2): the resulting EL_N had a genuine addot
dependence the paper's be1 does not have, and directly checking the
claimed IBP identity's residual confirmed it is non-zero. The current
code computes the G_{mu nu} nabla^mu A^nu piece DIRECTLY from the actual
Einstein tensor (via tensor_utils.compute_curvature(), the same
machinery fR_gravity.py's covariant cross-check route already uses) --
slower than either IBP shortcut, but verified: residual against
arXiv:1703.09573's be1 is exactly 0 for a G5-only truncated theory and
for the full random-multi-theory sweep in
sympy_layer/proca_background_1703_09573.py.

STATUS (run_checks()): with the direct-Einstein-tensor L5 G-piece AND
the corrected Euler-Lagrange formula for N (see euler_lagrange_all()'s
docstring -- EL_N previously used the naive dL/dN, missing a genuine
d/dt(dL/dNdot) contribution), all checks in this file plus the H1
cross-check against arXiv:1703.09573's be1/be2/be3 in
sympy_layer/proca_background_1703_09573.py pass. The former "KNOWN OPEN
ISSUE" (residual addot in the A0 equation with L5 included) is
RESOLVED -- its root cause was the flawed L5 G-piece reduction above,
not (as originally suspected) an error in how EL_A0 combines the
gravity-piece and cubic-bracket pieces.
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


def build_lagrangian(include_L5=True, d2_val=None,
                      G2_expr=None, G3_expr=None, G4_expr=None, G5_expr=None):
    """G2_expr..G5_expr: optional concrete sympy expressions in the module's
    Xsym=sp.Symbol('X', positive=True) (independently re-declared with the
    same name+assumptions is fine -- sympy interns Symbol('X',
    positive=True) so it compares equal/substitutable across call sites,
    verified before relying on this). Default None keeps the original
    abstract sp.Function('Gi')(Xsym) behavior (unchanged, so existing
    regression is preserved). Concrete forms let a caller (e.g.
    verification/paper_equations/arxiv_1703_09573.py's H1 comparison)
    numerically evaluate G_i,X / G_i,XX without fighting sympy's
    abstract-UndefinedFunction substitution machinery."""
    coords, g, ginv, Gamma, a_t, N_t = flrw_geometry()
    A0_t = sp.Function('A0')(t)
    A_lower = [A0_t, 0, 0, 0]

    nabla_A, M = covariant_derivative_A(Gamma, ginv, A_lower, coords)
    div_A = sp.simplify(sum(M[i, i] for i in range(4)))                 # nabla_mu A^mu
    MM_trace = sp.simplify(sum(M[i, j] * M[j, i] for i in range(4) for j in range(4)))  # nabla_mu A_nu nabla^nu A^mu

    # X = -1/2 A_mu A^mu = -1/2 g^{00} A0^2 = A0^2/(2N^2)
    X_t = sp.simplify(-sp.Rational(1, 2) * ginv[0, 0] * A0_t**2)

    Xsym = sp.Symbol('X', positive=True)
    G2 = G2_expr if G2_expr is not None else sp.Function('G2')(Xsym)
    G3 = G3_expr if G3_expr is not None else sp.Function('G3')(Xsym)
    G4 = G4_expr if G4_expr is not None else sp.Function('G4')(Xsym)
    G5 = G5_expr if G5_expr is not None else sp.Function('G5')(Xsym)
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
    # addot coefficient, integrate by parts once). That FAILED (the
    # extracted coefficient itself depends on adot, so differentiating it
    # re-introduces addot via the chain rule; did not converge after 10
    # passes -- a genuine dead end, not just slow convergence).
    #
    # A SECOND attempt used the contracted Bianchi identity
    # nabla_mu G^{mu nu} = 0 to rewrite this as a boundary term (dropped)
    # plus G5,X Xdot g^{00} A0 -- documented at length in an earlier
    # version of this function as "the working approach". IT WAS NOT: an
    # H1 comparison against arXiv:1703.09573's own be1 (see
    # verification/paper_equations/arxiv_1703_09573.py and
    # sympy_layer/proca_background_1703_09573.py) exposed that the
    # resulting EL_N had a genuine, non-vanishing addot dependence that
    # the paper's be1 does not have. Directly checking the claimed IBP
    # identity's residual (full symbolic expansion, not just the resulting
    # EL_N) confirmed it is NON-ZERO -- the "boundary term is pure
    # nabla_mu(...)" step in that derivation has an error (most likely:
    # the claim that V^mu=G5 G^{mu nu} A_nu is purely temporal, hence
    # sqrt(-g) nabla_mu V^mu is an EXACT total time derivative, does not
    # survive contact with the actual covariant divergence formula
    # sqrt(-g) nabla_mu V^mu = partial_mu(sqrt(-g) V^mu) once genuine N(t)
    # dependence inside the Einstein tensor itself is accounted for --
    # not re-derived by hand a third time here, since a direct computation
    # sidesteps the question entirely and was independently verified below).
    #
    # WORKING APPROACH (verified, not shortcut): compute the Einstein
    # tensor G_{mu nu} = R_{mu nu} - (1/2) g_{mu nu} R directly from the
    # SAME Ricci/Ricci-scalar machinery already used and cross-validated
    # elsewhere in this project (tensor_utils.compute_curvature(), the
    # same function fR_gravity.py's covariant route B relies on), for
    # the metric with a fully general, un-gauge-fixed N(t) -- exactly the
    # "具體分量計算只保留給最後的 ansatz 代入" (component calculation
    # reserved for the final ansatz substitution) principle docs/SPEC.md
    # module J calls for, since this IS the final FLRW+temporal-A ansatz
    # substitution, not a general covariant derivation. Verified against
    # arXiv:1703.09573's be1 exactly (H1; residual identically 0 for a
    # G5-only truncated theory, and the full random-multi-theory sweep in
    # sympy_layer/proca_background_1703_09573.py's run_checks()).
    from tensor_utils import compute_curvature
    _, _, Ric, Rs = compute_curvature(g, coords)
    Einstein_lower = sp.Matrix(n, n, lambda i, j: sp.simplify(Ric[i, j] - sp.Rational(1, 2) * g[i, j] * Rs))
    L5_Gpiece_density = sp.simplify(sum(Einstein_lower[mu, nu] * N_up2[mu, nu]
                                         for mu in range(n) for nu in range(n)))
    L5_Gpiece_reduced = sp.expand(N_t * a_t**3 * G5_t * L5_Gpiece_density)

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

    # EL w.r.t. N: constraint ("00" equation).
    #
    # BUG FIX (found via H1 comparison against arXiv:1703.09573's own be1/
    # be2 -- see verification/paper_equations/arxiv_1703_09573.py and
    # sympy_layer/proca_background_1703_09573.py): this was previously
    # `EL_N = sp.diff(L, N_t)`, i.e. the NAIVE dL/dN alone. That is only
    # correct if L has no explicit Ndot dependence. It does: L3_reduced,
    # L4_R_reduced (via Fdot = G4,X(X) Xdot) and L5_Gpiece_reduced (via
    # Xdot) all involve Xdot = d/dt(A0^2/(2N^2)), which contains Ndot
    # through the chain rule whenever N(t) is kept general (confirmed via
    # `L.has(sp.diff(N_t, t))` before this fix). Since N is a genuine
    # independent generalized coordinate of this reduced mini-superspace
    # Lagrangian, its correct Euler-Lagrange equation MUST use the same
    # general form used below for a and A0: d/dt(dL/dqdot) - dL/dq = 0.
    # An earlier version of this fix used the opposite sign order
    # (dL/dN - d/dt(dL/dNdot)) and, while it "looked" plausible in
    # isolation, it disagreed with the a/A0 sign convention and produced
    # a G2-term-vs-G3-term RELATIVE sign inconsistency against the paper's
    # be1 that only showed up once compared term-by-term (a single
    # overall sign flip of EL_N can never fix a relative-sign
    # inconsistency between two different G_i's contributions -- both
    # orderings were tried and compared explicitly before settling on
    # this one). Verified: with this order AND the correct phi=-A0
    # index-convention mapping (see
    # verification/paper_equations/arxiv_1703_09573.py's module
    # docstring), EL_N = -a^3 * be1 exactly, uniformly across G2, G3, G4
    # terms tested independently (sympy_layer/proca_background_1703_09573.py).
    dL_dNdot = sp.diff(L, sp.diff(N_t, t))
    EL_N = sp.diff(dL_dNdot, t) - sp.diff(L, N_t)

    # EL w.r.t. a: dynamical ("ii") equation.
    #
    # BUG FIX (found the same way as the EL_N fix above, via H1 comparison
    # -- see docs/SPEC.md Stage 2): this was previously the plain
    # first-order `d/dt(dL/dadot) - dL/da`, which is only correct if L has
    # no explicit addot dependence. Since the L5 G-piece fix (see
    # build_L5_term's docstring) now uses the actual Einstein tensor, L
    # DOES depend on addot (confirmed via `L.has(sp.diff(a_t,t,2))`) --
    # the Ricci scalar/tensor genuinely contains second time-derivatives
    # of a(t). For a Lagrangian L(a, adot, addot), the correct
    # (Ostrogradsky) Euler-Lagrange equation is
    #   dL/da - d/dt(dL/dadot) + d^2/dt^2(dL/daddot) = 0
    # i.e., matching this file's sign convention (d/dt(dL/dq) - dL/dq = 0
    # form used for N and A0 above):
    #   EL_a = d/dt(dL/dadot) - dL/da - d^2/dt^2(dL/daddot)
    # Verified: with this correction, EL_a = -3 a^2 be2 exactly (same
    # constant as the L2+L3+L4-only, no-addot-dependence case), across
    # independent random G2..G5 (sympy_layer/proca_background_1703_09573.py).
    # Missing this term previously produced spurious THIRD time-derivative
    # terms in EL_a (Derivative(a(t),(t,3))) that have no business
    # appearing in a background equation for a theory constructed
    # specifically to keep the equations of motion second-order -- a
    # useful sanity signal that something was missing, in hindsight.
    dL_dadot = sp.diff(L, sp.diff(a_t, t))
    dL_da = sp.diff(L, a_t)
    dL_daddot = sp.diff(L, sp.diff(a_t, t, 2))
    EL_a = sp.diff(dL_dadot, t) - dL_da - sp.diff(dL_daddot, t, 2)

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

    # Check: L_total DOES now legitimately contain addot/Adotdot before EL
    # is taken -- EXPECTED and CORRECT now that the L5 G-piece is computed
    # directly from the real Einstein tensor (see build_L5_term's
    # docstring) rather than the old, flawed "addot-free by construction"
    # reduced shortcut. A raw mini-superspace Lagrangian is allowed to
    # contain addot; what must be addot-free is the RESULTING equations of
    # motion where the paper says so (checks 4 and 5 below) -- an earlier
    # version of this check asserted the opposite (L_total itself must be
    # addot-free) as a proxy for "the reduction worked", which was only
    # ever true because that reduction was itself buggy (see
    # docs/SPEC.md Stage 2 and this module's docstring).
    L = data['L_total']
    ok4 = L.has(addot) or L.has(Addot0)
    results.append(('3. [with L5] L_total legitimately contains addot/Adotdot pre-EL (expected: the '
                     'direct-Einstein-tensor L5 G-piece is not artificially addot-free; only the '
                     'resulting EL_N/EL_A0 need to be, per checks 4-5)',
                     ok4, None))

    # Check: EL_N (Hamiltonian/Friedmann-like constraint) with L5 included
    # is still addot-free -- matches the paper's eq (2.12)/(be1), which
    # has no time derivative of phi or a beyond H itself.
    ok5 = not EL_N.has(addot)
    results.append(('4. [with L5] EL_N (Friedmann-like eq) is addot-free, matching eq (2.12) structure',
                     ok5, EL_N))

    # ---- RESOLVED (was "KNOWN OPEN ISSUE" in an earlier version of this
    # file): EL_A0 with L5 included previously had a residual addot term,
    # contradicting the paper's eq (be3), which is COMPLETELY algebraic
    # (no time derivatives of phi OR a) for arbitrary G2-G5. Root cause
    # found via H1 comparison (docs/SPEC.md Stage 2): the L5 gravity-piece
    # reduction (the old "Bianchi identity" shortcut) was itself wrong,
    # not how EL_A0 combined it with the cubic-bracket piece as originally
    # suspected. Fixed in build_L5_term() (direct Einstein-tensor
    # computation); this check now passes.
    dep_on_addot_A0 = EL_A0.has(addot)
    results.append(('5. [with L5] A0 equation (paper eq be3) is addot-free -- RESOLVED, was a known '
                     'open issue caused by a bug in the L5 gravity-piece reduction, now fixed',
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

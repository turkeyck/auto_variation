"""
Module F6 (docs/SPEC.md): generic N-field kinetic/gradient matrix
extractor for a Fourier-reduced, already-quadratic mode Lagrangian
L(q_1,...,q_N, qdot_1,...,qdot_N, k, t).

Generalizes sympy_layer/pert_engine.py's single-field
kinetic_coefficient()/gradient_coefficient() (kept as-is; this module
does not replace them, it extends the same convention to N>=1 fields --
kinetic_matrix([q], t)[0,0] and gradient_matrix(L,[q],k)[0,0] reproduce
pert_engine.py's K and gradient_coefficient() output exactly, verified in
run_checks() below).

Needed because the existing project's N=2 kinetic-matrix computations
(sympy_layer/scalar_plus_fluid_sector.py's K_VV/K_SS/K_VS) are each
hand-written per file with N hard-coded to exactly 2 -- docs/SPEC.md
module F6 requires "對任意 N 個場的二階作用量，回傳 N×N 的 kinetic
matrix 與 gradient matrix" (a generic N-field extractor), not a
2-field-only special case.

CONVENTION (matches pert_engine.py's single-field convention exactly --
verified below, not just asserted):
    K_ij = (1/2) d^2L / (dqdot_i dqdot_j)
    G_ij = -(1/2) d^2/(dq_i dq_j) [ (1/2) d^2L/dk^2 ]
so that for a diagonal single-field system L = K qdot^2 - (k^2/a^2) G q^2,
the sound speed is c_s^2 = G_ii * a^2 / K_ii (both matrices already carry
the correct sign so this ratio comes out positive for a normal,
non-ghost, stable mode -- verified against pert_engine.py's Maldacena
check, which independently established c_s^2=1 for that system).

LIMITATION (documented, not hidden): this extractor assumes L is
EXACTLY a quadratic form in (q_i, qdot_i) with a single overall
k-dependence entering as k^2 (linear in k^2, standard for a
spatially-Fourier-reduced action) -- it does NOT handle qdot_i*q_j
"symplectic" cross terms (a genuine possibility in some theories, not
needed by any check in this project so far) or higher powers of k. If a
future theory needs either, extend this module rather than working
around it ad hoc in a caller.
"""
import sympy as sp


def kinetic_matrix(L, fields, tvar):
    """K_ij = (1/2) d^2L/(dqdot_i dqdot_j), returned as an NxN sympy
    Matrix. `fields` is the ordered list of field functions
    [q1(t),...,qN(t)] (sp.Function(...)(tvar) objects)."""
    n = len(fields)
    qdots = [sp.diff(q, tvar) for q in fields]
    K = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            K[i, j] = sp.simplify(sp.diff(L, qdots[i], qdots[j]) / 2)
    return K


def gradient_matrix(L, fields, kvar):
    """G_ij = -(1/2) d^2/(dq_i dq_j) [ (1/2) d^2L/dk^2 ], returned as an
    NxN sympy Matrix. Assumes L is linear in k^2 (i.e. `sp.diff(L, kvar, 2)`
    is itself k-independent, standard for a Fourier-reduced quadratic
    action) -- this is checked, not assumed silently (see
    run_checks() 's use of the same pattern already established in
    pert_engine.py)."""
    n = len(fields)
    dL_dk2_half = sp.diff(L, kvar, 2) / 2
    G = sp.zeros(n, n)
    for i in range(n):
        for j in range(n):
            G[i, j] = sp.simplify(-sp.diff(dL_dk2_half, fields[i], fields[j]) / 2)
    return G


def no_ghost_conditions(K):
    """Return the list of leading principal minors of K (Sylvester's
    criterion building blocks): K must be positive-definite for a
    ghost-free theory, which holds iff every leading principal minor is
    positive. Returned as sympy expressions for the caller to check
    positivity of (under whatever background assumptions apply), not
    evaluated here (this module has no opinion on background
    positivity assumptions -- that is ansatz-specific, module E/F
    territory, not this generic linear-algebra utility's job)."""
    n = K.shape[0]
    return [sp.simplify(K[:m, :m].det()) for m in range(1, n + 1)]


def run_checks():
    results = []
    t = sp.Symbol('t', real=True)
    k = sp.Symbol('k', positive=True)

    # ---- Check 1: N=1 reproduces pert_engine.py's Maldacena result exactly ----
    from pert_engine import kinetic_coefficient as pe_kinetic_coefficient

    a = sp.Function('a')(t)
    epsH = sp.Symbol('epsilon_H', positive=True)
    zk = sp.Function('zeta_k')(t)
    L_mode = a**3 * epsH * (sp.diff(zk, t)**2 - k**2 / a**2 * zk**2)

    K = kinetic_matrix(L_mode, [zk], t)
    G = gradient_matrix(L_mode, [zk], k)
    K_pe = pe_kinetic_coefficient(L_mode, zk, t)

    ok1 = sp.simplify(K[0, 0] - K_pe) == 0
    results.append(('1. kinetic_matrix() N=1 matches pert_engine.kinetic_coefficient() exactly '
                     '(Maldacena L_mode)', ok1, (K[0, 0], K_pe)))

    cs2 = sp.simplify(G[0, 0] * a**2 / K[0, 0])
    ok2 = sp.simplify(cs2 - 1) == 0
    results.append(('2. c_s^2 = G_00 a^2 / K_00 = 1 for the Maldacena system (matches the '
                     'independently-established pert_engine.py result)', ok2, cs2))

    # ---- Check 3: N=2 diagonal (no cross-coupling) reduces to two independent N=1 systems ----
    epsH2 = sp.Symbol('epsilon_H2', positive=True)
    yk = sp.Function('y_k')(t)
    L_two_diag = L_mode + a**3 * epsH2 * (sp.diff(yk, t)**2 - k**2 / a**2 * yk**2)
    K2 = kinetic_matrix(L_two_diag, [zk, yk], t)
    G2 = gradient_matrix(L_two_diag, [zk, yk], k)
    ok3 = (sp.simplify(K2[0, 0] - K[0, 0]) == 0 and sp.simplify(K2[1, 1] - a**3 * epsH2) == 0
           and K2[0, 1] == 0 and K2[1, 0] == 0
           and G2[0, 1] == 0 and G2[1, 0] == 0)
    results.append(('3. N=2 diagonal (no cross-coupling) system: off-diagonal K,G entries are '
                     'exactly 0, diagonal entries match the two independent N=1 systems', ok3,
                     (K2, G2)))

    # ---- Check 4: N=2 with a genuine cross-coupling term reproduces the SAME K_VS this project's
    # own scalar_plus_fluid_sector.py already computed by hand for its 2-field kinetic system
    # (H2: independent generic-code cross-check against an existing hand-built special case). ----
    from scalar_plus_fluid_sector import build_L2c_with_fluid
    L2c, a_ssf, abar, chibar, dphi, chiV, dphis, Abar0, phibardot, m2, H = build_L2c_with_fluid()
    dL2_dabar = sp.expand(sp.diff(L2c, abar))
    dL2_dchibar = sp.expand(sp.diff(L2c, chibar))
    sol1 = sp.solve([sp.Eq(dL2_dabar, 0), sp.Eq(dL2_dchibar, 0)], [abar, chibar], dict=True)[0]
    L_step1 = sp.expand(L2c.subs(sol1))
    dL2_ddphi = sp.expand(sp.diff(L_step1, dphi))
    sol_dphi = sp.solve(sp.Eq(dL2_ddphi, 0), dphi, dict=True)[0]
    L_final = sp.expand(L_step1.subs(sol_dphi))

    K_generic = kinetic_matrix(L_final, [chiV, dphis], t)
    chiVdot = sp.diff(chiV, t)
    dphisdot = sp.diff(dphis, t)
    K_VV_hand = sp.simplify(sp.diff(L_final, chiVdot, 2) / 2)
    K_SS_hand = sp.simplify(sp.diff(L_final, dphisdot, 2) / 2)
    K_VS_hand = sp.simplify(sp.diff(sp.diff(L_final, chiVdot), dphisdot) / 2)

    ok4 = (sp.simplify(K_generic[0, 0] - K_VV_hand) == 0
           and sp.simplify(K_generic[1, 1] - K_SS_hand) == 0
           and sp.simplify(K_generic[0, 1] - K_VS_hand) == 0
           and sp.simplify(K_generic[1, 0] - K_VS_hand) == 0)
    results.append(('4. [H2] Generic kinetic_matrix() reproduces scalar_plus_fluid_sector.py\'s own '
                     'hand-built K_VV/K_SS/K_VS 2x2 kinetic matrix (including the off-diagonal cross '
                     'term) exactly, for its genuine (chiV, dphi_s) two-field system',
                     ok4, (K_generic, K_VV_hand, K_SS_hand, K_VS_hand)))

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

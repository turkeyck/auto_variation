"""
Ground-truth check against arXiv:1703.09573v2 (De Felice, Heisenberg,
Tsujikawa, "Observational constraints on generalized Proca theories",
Phys. Rev. D 95, 123540).

Unlike every other check in this repo -- which tests INTERNAL consistency
(GR limits, no-ghost positivity, numeric spot-checks) -- this file
transcribes the paper's own equations verbatim and asks whether the code
reproduces them. That distinction matters: all 40 of the repo's internal
checks passed while two results were wrong. One of those (the lapse
equation, B3/B4) is now fixed in proca_minisuperspace.py; the other (the
off-shell scalar-sector background, D2) still stands.

Target equations (the ones the project is supposed to reproduce):
    (3.3)        tensor quadratic action,       via (3.4) qT, (3.5) cT^2
    (3.8)        vector sound speed,            via (3.7) qV
    (3.10)-(3.15) scalar-sector equations of motion, via (3.16)-(3.23) w1..w7
Supporting/background equations used to validate the foundations:
    (2.11)-(2.13) background equations of motion

CONVENTION (established by test A and confirmed independently by test B):
the repo's A0 is the LOWER-index component A_0, and the paper's phi is
the UPPER-index A^0, so
                          phi = -A0 .
Both give X = A0^2/2 = phi^2/2, so X-only expressions are insensitive to
this, but every term odd in phi (all the G3 terms, the G5,X terms in qT
and cT^2) flips sign under it. Nothing in the repo records this anywhere,
which is exactly how such a mismatch goes unnoticed.
"""
import sympy as sp

t = sp.Symbol('t', real=True)
k = sp.Symbol('k', positive=True)


# =====================================================================
# A.  paper (3.4) qT, (3.5) cT^2      vs   proca_tensor.py
# =====================================================================
def check_tensor_sector():
    import proca_tensor as PT

    phi, phidot, H = sp.symbols('phi phidot H', real=True)
    G4, G4X, G5X = sp.symbols('G4 G4X G5X', real=True)

    # --- paper, verbatim ---
    qT_paper = 2 * G4 - 2 * phi**2 * G4X + H * phi**3 * G5X          # (3.4)
    FT_paper = 2 * G4 + phi**2 * phidot * G5X                        # numerator of (3.5)

    # --- repo, with the phi = -A0 convention and X = phi^2/2 ---
    A0, A0dot, X = -phi, -phidot, phi**2 / 2
    qT_repo = PT.G_T(G4, G4X, X, A0, H, G5X)
    FT_repo = PT.F_T(G4, X, G5X, A0dot)

    out = []
    out.append(('A1. (3.4) qT matches proca_tensor.G_T  [under phi = -A0]',
                sp.simplify(qT_repo - qT_paper) == 0,
                sp.simplify(qT_repo - qT_paper)))
    out.append(('A2. (3.5) cT^2 numerator matches proca_tensor.F_T  [under phi = -A0]',
                sp.simplify(FT_repo - FT_paper) == 0,
                sp.simplify(FT_repo - FT_paper)))
    # the same formulas under the OTHER convention, to show the sign hazard is real
    qT_wrong = PT.G_T(G4, G4X, X, phi, H, G5X)
    out.append(('A3. HAZARD: feeding phi (not A0) into G_T gives the WRONG sign on the '
                'G5,X term -- the convention is undocumented in proca_tensor.py',
                sp.simplify(qT_wrong - qT_paper) != 0,
                sp.simplify(qT_wrong - qT_paper)))
    out.append(('A4. GAP: no code anywhere builds the tensor quadratic action (3.3) itself; '
                'qT/cT^2 are hard-coded coefficient formulas, not derived from the action',
                False, 'not implemented'))
    return out


# =====================================================================
# B.  paper (2.11)-(2.13)             vs   proca_minisuperspace.py
# =====================================================================
def _to_H(expr, a_t, H, Hd):
    e = expr.subs(sp.Derivative(a_t, (t, 2)), a_t * (H**2 + Hd))
    e = e.subs(sp.Derivative(a_t, t), a_t * H)
    return sp.expand(e)


def _proportional(repo_eq, paper_eq, forbidden):
    """repo_eq and paper_eq describe the same equation iff their ratio is a
    pure power of a(t) times a number (no H, no A0, no G-functions left)."""
    ratio = sp.simplify(sp.expand(repo_eq) / sp.expand(paper_eq))
    return (not any(ratio.has(s) for s in forbidden)), ratio


def check_background():
    import proca_minisuperspace as PM

    data = PM.build_lagrangian()
    L = data['L_total']
    a_t, N_t, A0_t = data['a_t'], data['N_t'], data['A0_t']
    Xsym, G2f, G3f, G4f = data['Xsym'], data['G2'], data['G3'], data['G4']

    H, Hd = sp.symbols('H Hdot', real=True)
    gauge = {N_t: 1, sp.diff(N_t, t): 0, sp.diff(N_t, t, 2): 0}
    Ndot = sp.diff(N_t, t)

    # ---- the repo's own EL equations, plus the naive one that drops -d/dt(dL/dNdot) ----
    EL_N_repo, EL_a, EL_A0 = PM.euler_lagrange_all(data)
    EL_N_naive = sp.simplify(sp.diff(L, N_t).subs(gauge))

    # ---- paper (2.11)-(2.13), verbatim, with G5 = 0 (repo has no L5) and rho_M = P_M = 0 ----
    X_t = data['X_t'].subs({N_t: 1})
    G2 = G2f.subs(Xsym, X_t)
    G2X = sp.diff(G2f, Xsym).subs(Xsym, X_t)
    G3X = sp.diff(G3f, Xsym).subs(Xsym, X_t)
    G4 = G4f.subs(Xsym, X_t)
    G4X = sp.diff(G4f, Xsym).subs(Xsym, X_t)
    G4XX = sp.diff(G4f, Xsym, 2).subs(Xsym, X_t)

    ph, phd = -A0_t, -sp.diff(A0_t, t)                 # phi = -A0

    eq211 = (G2 - G2X * ph**2 - 3 * G3X * H * ph**3 + 6 * G4 * H**2
             - 6 * (2 * G4X + G4XX * ph**2) * H**2 * ph**2)
    eq212 = (G2 - phd * ph**2 * G3X + 2 * G4 * (3 * H**2 + 2 * Hd)
             - 2 * G4X * ph * (3 * H**2 * ph + 2 * H * phd + 2 * Hd * ph)
             - 4 * G4XX * H * phd * ph**3)
    eq213 = ph * (G2X + 3 * G3X * H * ph + 6 * G4X * H**2 + 6 * G4XX * H**2 * ph**2)

    forbidden = [H, Hd, A0_t, G2f, G3f, G4f]
    out = []

    ok, r = _proportional(_to_H(EL_a, a_t, H, Hd), eq212, forbidden)
    out.append(('B1. (2.12) "ii" equation matches EL_a', ok, f'ratio = {r}'))

    ok, r = _proportional(_to_H(EL_A0, a_t, H, Hd), eq213, forbidden)
    out.append(('B2. (2.13) vector-field equation matches EL_A0', ok, f'ratio = {r}'))

    ok_repo, r_repo = _proportional(_to_H(EL_N_repo, a_t, H, Hd), eq211, forbidden)
    out.append(('B3. (2.11) Friedmann "00" equation matches EL_N', ok_repo, f'ratio = {r_repo}'))

    # Regression guard for the defect B3 used to expose: L_total genuinely depends on
    # Ndot (through the integrated-by-parts L3 term), so an EL_N built from dL/dN alone
    # gets the G3 sector wrong. Both halves must hold, or the fix has been reverted.
    ok_naive, r_naive = _proportional(_to_H(EL_N_naive, a_t, H, Hd), eq211, forbidden)
    out.append(('B4. REGRESSION GUARD: L depends on Ndot, and dropping -d/dt(dL/dNdot) '
                'would break (2.11) again',
                sp.simplify(sp.diff(L, Ndot)) != 0 and not ok_naive,
                f'dL/dNdot = {sp.simplify(sp.diff(L, Ndot))}; naive-EL_N ratio = {r_naive}'))

    out.append(('B5. GAP: L5 is absent from build_lagrangian(), so (2.11)-(2.13) can only be '
                'tested with G5 = 0 -- yet every G5 term in (3.4),(3.5),(3.8),(3.16)-(3.22) '
                'depends on it', False, 'not implemented'))
    return out


# =====================================================================
# C.  paper (3.7) qV, (3.8) cV^2      vs   vector_sector.py
# =====================================================================
def check_vector_sector():
    phi, phidot, H, Hdot, Mpl = sp.symbols('phi phidot H Hdot M_pl', positive=True)

    # the theory vector_sector.py actually implements: G2 = m^2 X + F, G4 = Mpl^2/2,
    # G3 = G5 = G6 = g5 = 0  ->  G2_F = 1 and every other coefficient vanishes
    G2F, G2Y, g5, G6, G6X, G4X, G5X = 1, 0, 0, 0, 0, 0, 0
    qT = Mpl**2

    qV = G2F + 2 * G2Y * phi**2 - 4 * g5 * H * phi + 2 * G6 * H**2 + 2 * G6X * H**2 * phi**2   # (3.7)
    cV2 = (1 + phi**2 * (2 * G4X - G5X * H * phi)**2 / (2 * qT * qV)
           + 2 * (G6 * Hdot - G2Y * phi**2
                  - (H * phi - phidot) * (H * phi * G6X - g5)) / qV)                           # (3.8)

    return [('C1. (3.8) cV^2 = 1 agrees with vector_sector.py -- but VACUOUSLY: every '
             'G-function-dependent term in (3.8) is identically zero for the implemented '
             'theory, so no structure of (3.8) is tested',
             False, f'specialised (3.7) qV = {sp.simplify(qV)}, (3.8) cV^2 = {sp.simplify(cV2)}')]


# =====================================================================
# D.  paper (3.10)-(3.15)             vs   scalar_sector.py
# =====================================================================
def check_scalar_sector():
    import scalar_sector as SS

    phi, H, Mpl, m2 = sp.symbols('phi H M_pl m2', positive=True)
    G4, qT, qV = Mpl**2 / 2, Mpl**2, 1

    # (3.16)-(3.22) and (3.26), specialised to GR + standard Proca
    # (G3 = G4,X = G5 = G6 = g5 = 0, G2 = m^2 X + F so G2,XX = 0)
    w1 = -4 * H * G4
    w2 = w1 + 2 * H * qT
    w4 = -6 * H**2 * G4
    w5 = w4 - sp.Rational(3, 2) * H * (w1 + w2)
    w6 = w7 = sp.Integer(0)
    qS = 3 * w1**2 + 4 * qT * w4

    degenerate = all(sp.simplify(x) == 0 for x in (w2, w5, w6, w7, qS))

    out = []
    out.append(('D1. (3.24) Q_S / (3.25) c_S^2 are usable for the implemented theory',
                not degenerate,
                f'w2={sp.simplify(w2)}, w5={sp.simplify(w5)}, w6={w6}, w7={w7}, '
                f'qS={sp.simplify(qS)} -> Q_S = 0 and c_S^2 = mu_S/(...qS) diverges'))

    # (2.13) for this theory reads phi * G2,X = m^2 phi = 0, forcing phi = 0.
    # (Textbook standard-Proca constraint: F^{mu 0} = 0 for a homogeneous A_0, so
    #  nabla_mu F^{mu nu} = m^2 A^nu gives m^2 A^0 = 0.)
    out.append(('D2. scalar_sector.py runs on an on-shell background',
                False,
                'background eq (2.13) forces phi = 0 for standard Proca, but build_L2c() '
                'carries a free non-zero Abar0(t): the background violates its own field equation'))

    out.append(('D3. (3.10)-(3.15) can be evaluated for the implemented theory',
                False,
                '(3.10),(3.12),(3.13),(3.23) all divide by phi, and phi = 0 is forced by D2 -- '
                'the paper\'s scalar parametrisation is singular at this operating point. '
                'Testing (3.10)-(3.15) requires switching on at least one of G3, G4(X), G5.'))

    # does the headline result survive once the background is made consistent?
    L2c, a, abar, chibar, dphi, chiV, Abar0, m2s, _ = SS.build_L2c()
    L0 = sp.expand(L2c.subs({Abar0: 0, sp.Derivative(Abar0, t): 0}))
    sol = sp.solve([sp.Eq(sp.diff(L0, abar), 0), sp.Eq(sp.diff(L0, chibar), 0)],
                   [abar, chibar], dict=True)
    L1 = sp.expand(L0.subs({abar: sol[0][abar], chibar: sol[0][chibar]}))
    sol_d = sp.solve(sp.Eq(sp.diff(L1, dphi), 0), dphi, dict=True)
    K = sp.simplify(sp.diff(L1.subs(dphi, sol_d[0][dphi]), sp.Derivative(chiV, t), 2) / 2)
    K_claimed = k**2 * m2s * a**3 / (2 * (k**2 + m2s * a**2))
    out.append(('D4. the reported K = k^2 m^2 a^3 / (2(k^2 + m^2 a^2)) survives fixing the '
                'background to Abar0 = 0 (so that result itself stands, even though the '
                'constraint system it was derived from should have been trivial)',
                sp.simplify(K - K_claimed) == 0,
                f'K(Abar0=0) = {K}; constraints give {sol[0]}'))
    return out


SECTIONS = [
    ('A. tensor sector -- paper (3.3)/(3.4)/(3.5)', check_tensor_sector),
    ('B. background    -- paper (2.11)/(2.12)/(2.13)', check_background),
    ('C. vector sector -- paper (3.7)/(3.8)', check_vector_sector),
    ('D. scalar sector -- paper (3.10)-(3.15)', check_scalar_sector),
]


if __name__ == '__main__':
    n_pass = n_total = 0
    for title, fn in SECTIONS:
        print(f'\n===== {title} =====')
        for name, ok, detail in fn():
            n_total += 1
            n_pass += bool(ok)
            print(f'[{"PASS" if ok else "FAIL"}] {name}')
            if not ok:
                print(f'       {detail}')
    print(f'\n{n_pass}/{n_total} checks passed')
    print("""
VERDICT: none of the target equations (3.3), (3.8), (3.10)-(3.15) is
reproduced by this repo.
  (3.3)        -- not built at all; only its coefficients qT, cT^2 are
                  hard-coded, and those DO match the paper (under phi = -A0).
  (3.8)        -- agrees only in the Maxwell limit, where every term that
                  makes (3.8) non-trivial has been set to zero.
  (3.10)-(3.15)-- cannot be compared: the implemented theory forces phi = 0,
                  where the paper's scalar variables are singular.
Two genuine defects were found that the repo's 40 internal checks all missed:
  * EL_N omitted -d/dt(dL/dNdot), corrupting the Friedmann equation when
    G3 != 0 -- FIXED in proca_minisuperspace.py; B3 now passes and B4 guards it
  * the scalar sector runs on a background that violates its own field
    equation (D2) -- STILL OPEN
""")

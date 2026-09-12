"""
Reproduce arXiv:1603.05806 equations (3.20)-(3.23) (labels LVM, Wiex,
the unlabeled u_i definition, and conser): the Schutz-Sorkin perfect-
fluid VECTOR-type perturbation sector.

CORRECTION: an earlier version of this file called these "(3.35)-
(3.38)" based on an ar5iv HTML auto-numbering that turned out to be
wrong. The actual numbers were confirmed by installing TeX Live and
running pdflatex on the paper's own arXiv e-print source, then reading
the equation numbers straight out of the resulting .aux file's
\\newlabel entries -- the authoritative source, not OCR, not an HTML
converter's guess, not memory. The physics/derivation below was never
affected by this -- only the citation label was wrong.

Source verified directly against the paper's own LaTeX (jcapresub1.tex,
downloaded from the arXiv e-print, not the OCR'd PDF -- exact, not
reconstructed from memory):

  Schutz-Sorkin action (paper eq. Spf):
    S_M = -int d^4x [ sqrt(-g) rho_M(n) + J^mu(d_mu ell + A_i d_mu B_i) ]
    n = sqrt(J^a J^b g_ab / g)

  Background+perturbation split (paper, vector sector):
    ell = -int rho_{M,n} dt'                    (v=0: no scalar velocity potential)
    J^0 = N0  (= n0 a^3, constant)               (delta J=0: no scalar density pert.)
    J^i = (1/a^2) W_i                            (delta j=0: only transverse W_i)
    A_i = delta A_i,  B_i = x_i + delta B_i       (Lin-constraint vector potentials)
    metric (vector part only): g_00=-1, g_0i=V_i, g_ij=a^2 delta_ij

This file builds S_M to quadratic order in (V_i, W_i, delta A_i, delta B_i)
directly from the Schutz-Sorkin action (not by copying the paper's
already-reduced result), and checks the result against paper eqs.
(3.20)-(3.23) term by term.
"""
import sympy as sp

t = sp.Symbol('t', real=True)
# No time derivatives of a(t) appear anywhere in this file (unlike the ADM
# files), so a plain positive symbol is used instead of Function(t) --
# this also lets sympy auto-simplify sqrt(a**6) -> a**3 (it can't do that
# for a generic Function(t), whose sign it doesn't know, which otherwise
# leaves spurious sqrt(a(t)**6) terms cluttering every comparison below).
a = sp.Symbol('a', positive=True)
rho_M = sp.Symbol('rho_M', real=True)
rho_Mn = sp.Symbol('rho_Mn', real=True)   # rho_{M,n} = d(rho_M)/dn, background value
N0 = sp.Symbol('N0', positive=True)  # curly-N_0 = n0 a^3, a constant


def build_action_from_scratch():
    """Build (S_V^(2))_M directly from the Schutz-Sorkin action, i.e.
    derive n to O(eps^2) from its own definition and expand rho_M(n)
    and the J.dell, J.A.dB pieces -- not copying eq (3.20).

    Uses explicit eps bookkeeping (every perturbation tagged eps*(...))
    and manual Taylor-coefficient extraction, NOT sequential single-
    variable sp.series() calls -- the latter fails to correctly discard
    cross terms like V_i^2 W_i^2 (each individually below the per-
    variable truncation degree, but jointly higher than quadratic
    order), which silently corrupted an earlier version of this file.
    """
    eps = sp.Symbol('eps_bk')
    Vi0, Wi0, dAi0 = sp.symbols('V_i W_i deltaA_i', real=True)
    dBidot0 = sp.Symbol('deltaBidot', real=True)  # d(delta B_i)/dt
    Vi, Wi, dAi, dBidot = eps * Vi0, eps * Wi0, eps * dAi0, eps * dBidot0

    # J^alpha J^beta g_alpha_beta = (J^0)^2 g_00 + 2 J^0 J^i g_0i + J^i J^i g_ii
    J0, Ji = N0, Wi / a**2
    JJg = J0**2 * (-1) + 2 * J0 * Ji * Vi + Ji**2 * a**2

    # g = det(g_munu) for g_00=-1, g_0i=V_i, g_ij=a^2 delta_ij (vector-only):
    # g = det(g_ij) * (g_00 - g_0i (g^ij) g_j0) = a^6 * (-1 - Vi^2/a^2)
    g_det = a**6 * (-1 - Vi**2 / a**2)

    n2 = sp.together(JJg / g_det)
    n0_val = N0 / a**3
    delta_n2 = sp.expand(sp.diff(n2, eps, 2).subs(eps, 0) / 2) * eps**2  # keep only the O(eps^2) piece
    # n = n0 sqrt(1 + delta_n2/n0^2) ~ n0 (1 + delta_n2/(2 n0^2)) to the
    # order needed (delta_n2 itself starts at O(eps^2), so this linear
    # sqrt expansion is exact to O(eps^2)).
    delta_n = sp.expand(delta_n2 / (2 * n0_val))

    sqrt_mg_full = sp.sqrt(-g_det)
    sqrt_mg = (a**3 + sp.diff(sqrt_mg_full, eps, 2).subs(eps, 0) / 2 * eps**2)

    # sqrt(-g) rho_M(n) ~ sqrt(-g)[rho_M(n0) + rho_Mn * delta_n]  (delta_n is O(eps^2))
    term_rho = sp.expand(sqrt_mg * (rho_M + rho_Mn * delta_n))
    term_rho_quad = sp.expand(sp.diff(term_rho, eps, 2).subs(eps, 0) / 2)

    # J^mu d_mu(ell): ell = -int rho_Mn dt' (background only for pure vector)
    # J^mu d_mu ell = J^0 * (-rho_Mn) = -N0 rho_Mn  (a pure background term, no eps dependence)

    # J^mu A_i d_mu(B_i) = dAi * [J^mu d_mu(x_i + dBi)] = dAi*[J^i (from d_mu x_i) + J^0 dBidot (leading order)]
    term_AB_full = dAi * (Ji + N0 * dBidot)
    term_AB_quad = sp.expand(sp.diff(term_AB_full, eps, 2).subs(eps, 0) / 2)

    L_full = -(term_rho_quad + term_AB_quad)
    L_full = sp.expand(L_full)
    return L_full, Vi0, Wi0, dAi0, dBidot0


def run_checks():
    results = []
    L_derived, Vi, Wi, dAi, dBidot = build_action_from_scratch()

    # ---- paper's eq (3.20), transcribed directly from its own bracket
    # (1/(2a^2 N0)){rho_Mn(W_i^2+N0^2 V_i^2) + N0(2 rho_Mn V_i W_i - a^3 rho_M V_i^2)}
    # - N0 deltaA_i deltaBidot - (1/a^2) W_i deltaA_i
    L_paper = sp.expand(
        rho_Mn * Wi**2 / (2 * a**2 * N0) + N0 * rho_Mn * Vi**2 / (2 * a**2)
        + rho_Mn * Vi * Wi / a**2 - a * rho_M / 2 * Vi**2
        - N0 * dAi * dBidot - Wi * dAi / a**2
    )

    diff = sp.simplify(L_derived - L_paper)
    ok1 = diff == 0
    results.append(('1. (S_V^(2))_M built from scratch (Schutz-Sorkin action, n from its own '
                     'definition) matches paper eq (3.20) term by term', ok1, diff))

    # ---- eq (3.21): vary w.r.t. W_i
    dL_dWi = sp.diff(L_derived, Wi)
    Wi_sol = sp.solve(sp.Eq(dL_dWi, 0), Wi)[0]
    Wi_sol = sp.simplify(Wi_sol)
    Wi_paper = N0 * (dAi - rho_Mn * Vi) / rho_Mn
    ok2 = sp.simplify(Wi_sol - Wi_paper) == 0
    results.append(('2. Varying w.r.t. W_i reproduces eq (3.21): W_i = N0(deltaA_i - rho_Mn V_i)/rho_Mn',
                     ok2, Wi_sol))

    # ---- eq (3.22)/(3.23): substitute W_i, vary w.r.t. delta A_i
    L_step1 = sp.expand(L_derived.subs(Wi, Wi_sol))
    dL_ddAi = sp.diff(L_step1, dAi)
    dAi_sol = sp.solve(sp.Eq(dL_ddAi, 0), dAi)[0]
    dAi_sol = sp.simplify(dAi_sol)
    # paper: delta A_i = rho_Mn * u_i,  u_i = V_i - a^2 * d(delta B_i)/dt
    ui_expected = Vi - a**2 * dBidot
    dAi_paper = rho_Mn * ui_expected
    ok3 = sp.simplify(dAi_sol - dAi_paper) == 0
    results.append(('3. Varying w.r.t. delta A_i (after substituting W_i) reproduces eq (3.22)+(delAi): '
                     'delta A_i = rho_Mn (V_i - a^2 deltaB_i_dot) = rho_Mn u_i', ok3, dAi_sol))

    # ---- eq (3.23): the conservation relation rho_Mn u_i = (rho_M+P_M)/n0 u_i = const
    # comes from varying w.r.t delta B_i (whose own EOM is the EL/time-derivative
    # equation for the now-substituted action -- check it reduces to d/dt(rho_Mn u_i)=0
    # i.e. rho_Mn u_i is a constant of motion).
    L_step2 = sp.expand(L_step1.subs(dAi, dAi_sol))
    # dBi doesn't appear in L_step2 except through dBidot (since only the
    # derivative entered); the EL equation for delta B_i is
    # -d/dt(dL/d(dBidot)) = 0 (no explicit dBi dependence) -> dL/d(dBidot) = const
    dL_ddBidot = sp.simplify(sp.diff(L_step2, dBidot))
    # dL/d(deltaBidot) = -N0 * rho_Mn * u_i; N0 is a genuine time-independent
    # constant (paper: "N0 = n0 a^3 is a constant"), so this being conserved
    # is equivalent to the paper's stated rho_Mn*u_i = const (eq 3.23).
    ok4 = sp.simplify(dL_ddBidot - (-N0 * rho_Mn * ui_expected)) == 0
    results.append(('4. EL derivative w.r.t. delta B_i involves only d(dL/d(deltaBidot))/dt=0, i.e. '
                     'dL/d(deltaBidot) proportional to rho_Mn*u_i is conserved, matching eq (3.23)',
                     ok4, dL_ddBidot))

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

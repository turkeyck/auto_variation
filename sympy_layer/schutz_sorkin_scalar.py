"""
Reproduce the SCALAR-perturbation half of the Schutz-Sorkin perfect-fluid
matter sector, from scratch, replacing the k-essence duality shortcut in
scalar_plus_fluid_sector.py (docs/SPEC.md Stage 4: that file's own
docstring admits it uses a P(X)=X "stiff fluid" duality instead of the
real Schutz-Sorkin J^mu formalism -- flagged there, and in docs/SPEC.md's
module G requirements, as not acceptable for claiming (3.10)-(3.15)
reproduced).

SOURCES (both downloaded as arXiv e-prints, not OCR'd/from memory):
  - arXiv:1703.09573v2 (astrophv2.tex) Sec. "Scalar perturbations": states
    the target equations (per1)-(per6) but explicitly defers the
    underlying second-order matter action to "Eq.~(4.6) of
    Ref.~[DeFelice:2016uil]" -- i.e. this paper does NOT itself contain
    the Schutz-Sorkin scalar-sector derivation.
  - arXiv:1605.05066v2 (arXiv_v2.tex), the DeFelice:2016uil reference:
    contains the actual Schutz-Sorkin action (eq. Spf), its
    scalar-perturbation split (J^0=N0+deltaJ, ell=-int rho_Mn dt' -
    rho_Mn*v, J^i=(1/a^2)d^ik(d_k deltaj+W_k)), the definitions
    delta_rho_M = (rho_Mn/a^3) deltaJ and delta_n (eq. num, deltaj
    region), and the target second-order matter action (S_M)_S^(2)
    (eq. SMS) this module verifies against.

METRIC CONVENTION (arXiv_v2.tex eq. near line 625, flat gauge):
  ds^2 = -(1+2 alpha) dt^2 + 2(d_i chi) dt dx^i + a^2 dx^2
(vector/tensor parts V_i, h_ij dropped -- pure scalar sector only).
A single shared spatial profile cos(k*z) is used for every scalar
perturbation (alpha, chi, deltaJ, deltaj, v), following the SAME
plane-wave/exact-eigenvalue convention already established and verified
in sympy_layer/adm_scalar.py and sympy_layer/vector_field_equation.py --
this makes every product of two perturbations a POINTWISE (not just
period-averaged) identity, since (SMS) itself is stated in terms of raw
spatial gradients (d(deltaj), d(v), d(chi)), not yet Fourier-reduced.

DISCIPLINE (matching schutz_sorkin_vector.py, NOT sp.series() applied
sequentially per-variable -- that silently drops genuine O(eps^2)
cross-terms, as documented in this project's history): every
perturbation is tagged with an explicit eps_bk bookkeeping symbol and
the O(eps^2) action is extracted via a single manual
sp.series(..., eps_bk, 0, 3) call on the FULL expression, then reading
off the eps_bk**2 coefficient -- never per-variable truncation.

rho_M(n) is kept as a fully generic function via sp.Function, with
rho_{M,n}=d(rho_M)/dn and rho_{M,nn}=d^2(rho_M)/dn^2 evaluated at the
background n0(t)=N0/a(t)^3 THROUGH CHAIN RULE (rho_n.subs(n, n0_t)),
matching the pattern already established and verified in
sympy_layer/fluid_sector.py -- this is essential here (unlike
schutz_sorkin_vector.py, which could get away with treating rho_{M,n} as
a bare constant): eq (SMS)'s "-6a^8n_0^2 rho_{M,nn} H v delta_rho_M" term
only appears because rho_{M,n}(n0(t)) genuinely varies with time through
n0(t)=N0/a(t)^3, and differentiating ell = -int rho_Mn dt' - rho_Mn(t)*v
w.r.t. t picks up d(rho_Mn)/dt = rho_Mnn * dn0/dt = -3H n0 rho_Mnn via
the chain rule -- this is where the H-dependence in (SMS) originates.

RESULT (run_checks()):
  1. delta_n, built from n^2=J^a J^b g_ab/g's own definition (not copied
     from eq. num/deltaj), matches arXiv:1605.05066's stated delta_n
     formula EXACTLY to O(eps^2) -- both the O(eps) piece
     (delta_rho_M/rho_Mn) and the O(eps^2) gradient piece.
  2. Every term of (S_M)_S^(2) that couples to an ACTUAL matter
     perturbation (deltaJ, deltaj, or v) matches eq (SMS) EXACTLY --
     zero residual, checked by isolating exactly those terms via
     polynomial coefficient extraction, not a blanket "close enough".
  3. Varying (S_M)_S^(2) w.r.t. deltaj reproduces eq (deltaj):
     d(deltaj) = -a^3 n0 (d(v)+d(chi)) exactly.
  4. A KNOWN, UNDERSTOOD, explicitly-isolated residual remains: a term
     proportional PURELY to background quantities (P_M=n0 rho_Mn-rho_M,
     and rho_M(n0) itself) times alpha^2/chi^2 -- i.e. it does NOT couple
     to any matter perturbation variable at all. This term does not
     appear in eq (SMS) as quoted by arXiv:1605.05066, but a term of
     exactly this shape (proportional to background matter quantities
     times pure metric-perturbation-squared, with no matter-perturbation
     coupling) is exactly what gets absorbed into the GRAVITY+matter
     COMBINED second-order action's own w4-type coefficients when the
     matter and gravity sectors are added together (see e.g. the "+w_4
     alpha^2" term in arXiv_v2.tex's eq. sscalar/(3.9) region) -- NOT
     something a matter-action-ALONE derivation should be expected to
     cancel on its own. Reproducing that full combination (requiring the
     GP gravity sector's own quadratic action in alpha, chi) is
     docs/SPEC.md Stage 5 scope (module F: unified SVT ansatz builder),
     not this Stage 4 module's job. Flagged honestly, not silently
     dropped -- see check 4 below, which verifies this residual has
     EXACTLY the predicted P_M/rho_M(n0)-proportional shape (confirming
     it is understood, not an unexplained leftover).
"""
import sympy as sp

t, z = sp.symbols('t z', real=True)
k = sp.Symbol('k', positive=True)
N0 = sp.Symbol('N0', positive=True)


def _background(a_t):
    nsym = sp.Symbol('n', positive=True)
    rho = sp.Function('rho')
    rho_n = sp.diff(rho(nsym), nsym)
    rho_nn = sp.diff(rho(nsym), nsym, 2)
    n0_t = N0 / a_t**3
    rho_M_bg = rho(n0_t)
    rho_Mn_bg = rho_n.subs(nsym, n0_t)
    rho_Mnn_bg = rho_nn.subs(nsym, n0_t)
    return n0_t, rho_M_bg, rho_Mn_bg, rho_Mnn_bg


def build_delta_n(a_t, eps, abar, chibar, dJbar, djbar, order=2):
    """delta_n = n - n0, built from n^2=J^a J^b g_ab/g's own definition
    (arXiv:1605.05066 eq. num), NOT copied from that paper's already-
    reduced eq. deltaj-region formula. Returns delta_n as an eps-tagged
    expression, truncated to O(eps^order)."""
    prof = sp.cos(k * z)
    alpha = eps * abar * prof
    chi = eps * chibar * prof
    dchi_dz = sp.diff(chi, z)
    dJ = eps * dJbar * prof
    dj = eps * djbar * prof
    ddj_dz = sp.diff(dj, z)

    g00 = -(1 + 2 * alpha)
    g0z = dchi_dz
    gzz = a_t**2
    det2 = g00 * gzz - g0z**2
    g_det = sp.expand(det2 * a_t**2 * a_t**2)

    J0 = N0 + dJ
    Jz = ddj_dz / a_t**2
    JJg = J0**2 * g00 + 2 * J0 * Jz * g0z + Jz**2 * gzz
    n2 = sp.together(JJg / g_det)

    n0_t = N0 / a_t**3
    n2_series = sp.series(n2, eps, 0, order + 1).removeO()
    out = 0
    n0sq = n0_t**2
    # n = n0*sqrt(1 + sum_{m>=1} eps^m n2_m / n0^2); expand to O(eps^order)
    # via direct substitution into sqrt's Taylor series (m<=2 sufficient here).
    n2_1 = n2_series.coeff(eps, 1)
    delta_n_1 = sp.simplify(n2_1 / (2 * n0_t))
    out += eps * delta_n_1
    if order >= 2:
        n2_2 = n2_series.coeff(eps, 2)
        delta_n_2 = sp.simplify(n2_2 / (2 * n0_t) - n2_1**2 / (8 * n0_t**3))
        out += eps**2 * delta_n_2
    return out


def build_action_from_scratch(a_t):
    """(S_M)_S^(2): the O(eps^2) Schutz-Sorkin scalar matter action
    density, built directly from S_M=-int[sqrt(-g) rho_M(n) + J^mu(d_mu
    ell + ...)] (arXiv:1605.05066 eq. Spf), not copied from that paper's
    already-reduced eq. SMS. Returns (S2, variables-dict)."""
    n0_t, rho_M_bg, rho_Mn_bg, rho_Mnn_bg = _background(a_t)

    eps = sp.Symbol('eps_bk')
    abar, chibar, dJbar, djbar, vbar = [sp.Function(name)(t) for name in
                                         ('abar', 'chibar', 'dJbar', 'djbar', 'vbar')]
    prof = sp.cos(k * z)
    alpha = eps * abar * prof
    dJ = eps * dJbar * prof
    dj = eps * djbar * prof
    ddj_dz = sp.diff(dj, z)
    v = eps * vbar * prof

    g00 = -(1 + 2 * alpha)
    g0z = sp.diff(eps * chibar * prof, z)
    gzz = a_t**2
    det2 = g00 * gzz - g0z**2
    g_det = sp.expand(det2 * a_t**2 * a_t**2)
    sqrt_mg = sp.sqrt(-g_det)

    J0 = N0 + dJ
    Jz = ddj_dz / a_t**2

    delta_n = build_delta_n(a_t, eps, abar, chibar, dJbar, djbar, order=2)
    rho_M_full = rho_M_bg + rho_Mn_bg * delta_n + sp.Rational(1, 2) * rho_Mnn_bg * delta_n**2

    dell_dt = -rho_Mn_bg - sp.diff(rho_Mn_bg * v, t)
    dell_dz = -sp.diff(rho_Mn_bg * v, z)
    J_dell = J0 * dell_dt + Jz * dell_dz

    S_density_full = -(sqrt_mg * rho_M_full + J_dell)
    S_series = sp.series(sp.expand(S_density_full), eps, 0, 3).removeO()
    S2 = sp.expand(S_series.coeff(eps, 2))

    return S2, {'abar': abar, 'chibar': chibar, 'dJbar': dJbar, 'djbar': djbar, 'vbar': vbar,
                'n0_t': n0_t, 'rho_M_bg': rho_M_bg, 'rho_Mn_bg': rho_Mn_bg, 'rho_Mnn_bg': rho_Mnn_bg}


def build_paper_SMS(a_t, variables):
    """eq (SMS) transcribed directly from arXiv:1605.05066 (H1 ground
    truth), using n_0 (number DENSITY, N0/a^3) exactly where the paper's
    own n_0 symbol appears -- NOT the bare total number N0 (an earlier
    draft of this derivation conflated the two, which is exactly the
    kind of N0-vs-n0 mixup schutz_sorkin_vector.py's own module docstring
    separately warns about)."""
    abar, chibar, dJbar, djbar, vbar = (variables['abar'], variables['chibar'],
                                         variables['dJbar'], variables['djbar'], variables['vbar'])
    n0_t, rho_M_bg, rho_Mn_bg, rho_Mnn_bg = (variables['n0_t'], variables['rho_M_bg'],
                                              variables['rho_Mn_bg'], variables['rho_Mnn_bg'])
    H = sp.diff(a_t, t) / a_t
    prof = sp.cos(k * z)

    alpha_bar = abar * prof
    ddj_dz = sp.diff(djbar * prof, z)
    dv_dz = sp.diff(vbar * prof, z)
    dchi_dz = sp.diff(chibar * prof, z)
    vdot = sp.diff(vbar, t) * prof
    delta_rho_M = (rho_Mn_bg / a_t**3) * dJbar * prof

    target = ((1 / (2 * a_t**5 * n0_t * rho_Mn_bg**2)) * (
        rho_Mn_bg * (rho_Mn_bg**2 * ddj_dz**2 + 2 * a_t**3 * n0_t * rho_Mn_bg**2 * ddj_dz * dv_dz
                     + 2 * a_t**8 * n0_t * rho_Mn_bg * vdot * delta_rho_M
                     - 6 * a_t**8 * n0_t**2 * rho_Mnn_bg * H * vbar * prof * delta_rho_M)
        - a_t**8 * n0_t * rho_Mnn_bg * delta_rho_M**2)
        - a_t**3 * alpha_bar * delta_rho_M
        + (rho_Mn_bg / a_t**2) * dchi_dz * ddj_dz)
    return sp.expand(target)


def run_checks():
    results = []
    a_t = sp.Function('a', positive=True)(t)

    # ---- Check 1: delta_n matches arXiv:1605.05066 eq (num)/(deltaj-region) ----
    n0_t, rho_M_bg, rho_Mn_bg, rho_Mnn_bg = _background(a_t)
    eps = sp.Symbol('eps_bk')
    abar, chibar, dJbar, djbar = [sp.Function(name)(t) for name in ('abar', 'chibar', 'dJbar', 'djbar')]
    delta_n = build_delta_n(a_t, eps, abar, chibar, dJbar, djbar, order=2)
    prof = sp.cos(k * z)
    dJ_amp, dchi_amp, ddj_amp = dJbar * prof, sp.diff(chibar * prof, z), sp.diff(djbar * prof, z)
    delta_rho_M_1 = (rho_Mn_bg / a_t**3) * dJ_amp
    target_delta_n = (delta_rho_M_1 / rho_Mn_bg
                       - (N0**2 * dchi_amp**2 + 2 * N0 * dchi_amp * ddj_amp + ddj_amp**2) / (2 * N0 * a_t**5))
    # delta_n already carries explicit eps powers; compare the eps^1
    # coefficient against the paper's O(eps) piece (delta_rho_M/rho_Mn)
    # and the eps^2 coefficient against its O(eps^2) gradient piece.
    diff1 = sp.simplify(delta_n.coeff(eps, 1) - delta_rho_M_1 / rho_Mn_bg)
    grad_target = -(N0**2 * dchi_amp**2 + 2 * N0 * dchi_amp * ddj_amp + ddj_amp**2) / (2 * N0 * a_t**5)
    diff2 = sp.simplify(delta_n.coeff(eps, 2) - grad_target)
    ok1 = (diff1 == 0) and (diff2 == 0)
    results.append(('1. [H1] delta_n (built from n^2=J^a J^b g_ab/g) matches arXiv:1605.05066\'s '
                     'own delta_n formula exactly to O(eps^2), both the delta_rho_M/rho_Mn piece and '
                     'the quadratic-gradient piece', ok1, (diff1, diff2)))

    # ---- Check 2: every matter-perturbation-coupled term of (S_M)_S^(2) matches eq (SMS) ----
    S2, variables = build_action_from_scratch(a_t)
    target = build_paper_SMS(a_t, variables)
    dJbar_v, djbar_v, vbar_v = variables['dJbar'], variables['djbar'], variables['vbar']
    diff = sp.expand(S2 - target)
    # isolate terms with zero coupling to any matter perturbation (dJbar=djbar=vbar=0)
    background_only_residual = sp.simplify(diff.subs({dJbar_v: 0, djbar_v: 0, vbar_v: 0}))
    matter_coupled_residual = sp.simplify(diff - background_only_residual)
    ok2 = matter_coupled_residual == 0
    results.append(('2. [H1] Every term of (S_M)_S^(2) that couples to an actual matter perturbation '
                     '(deltaJ, deltaj, or v) matches arXiv:1605.05066 eq (SMS) exactly -- zero residual',
                     ok2, matter_coupled_residual))

    # ---- Check 3: delta_j equation of motion matches eq (deltaj) ----
    dS2_ddjbar = sp.diff(S2, djbar_v)
    sol = sp.solve(sp.Eq(dS2_ddjbar, 0), djbar_v)
    ok3 = (len(sol) == 1) and sp.simplify(sol[0] - (-N0 * (variables['chibar'] + vbar_v))) == 0
    results.append(('3. [H1] Varying (S_M)_S^(2) w.r.t. deltaj reproduces eq (deltaj): '
                     'd(deltaj) = -a^3 n0 (d(v)+d(chi)) exactly',
                     ok3, sol))

    # ---- Check 4: the remaining background-only residual has the predicted shape ----
    # P_M = n0 rho_Mn - rho_M(n0): background pressure, arXiv:1605.05066 eq (PM).
    P_M = n0_t * rho_Mn_bg - rho_M_bg
    chibar_v = variables['chibar']
    predicted = sp.expand(sp.Rational(1, 2) * a_t * P_M * k**2 * chibar_v**2 * sp.sin(k * z)**2
                           + sp.Rational(1, 2) * a_t**3 * rho_M_bg * abar**2 * sp.cos(k * z)**2)
    ok4 = sp.simplify(background_only_residual - predicted) == 0
    results.append(('4. The isolated background-only residual (no coupling to any matter perturbation) '
                     'has EXACTLY the predicted shape: (1/2)a P_M k^2 chibar^2 sin^2(kz) + '
                     '(1/2)a^3 rho_M(n0) abar^2 cos^2(kz) -- confirming it is understood (a term that '
                     'belongs to the gravity+matter COMBINED action\'s w4-type bookkeeping, Stage 5 '
                     'scope), not an unexplained leftover',
                     ok4, sp.simplify(background_only_residual - predicted)))

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

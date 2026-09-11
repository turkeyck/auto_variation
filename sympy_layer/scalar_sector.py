"""
Phase 3.4, steps 3-5 -- add Proca matter (mass term via X, Maxwell
kinetic term via F_{mu nu}F^{mu nu}) to the vacuum ADM scalar sector
from adm_scalar.py, solve the linear constraints for (abar, chibar) in
terms of the physical matter perturbations (dphi := delta A_0, chiV :=
longitudinal potential of delta A_i = partial_i chiV), and extract the
resulting kinetic structure. Target: check 3.e, GR + massive vector
field (standard Proca: Maxwell + mass term).

*** KNOWN ISSUE -- DO NOT TRUST THE RESULTS BELOW YET ***
Running this: the momentum constraint (dL2/dchibar = 0) comes out
depending ONLY on chiV (no chibar or abar self-coupling survives),
forcing chiV = 0 for generic parameters. That directly CONTRADICTS
the well-established physical fact that the Proca mass term is what
makes the vector field's longitudinal mode chiV genuinely dynamical
(the entire point of massive vs. massless vector cosmological
perturbation theory). So something in this file's construction is
wrong or incomplete -- most likely the hand-derived X_quad/F2_quad
expansions are missing a term, or there is a subtlety in how the flat
gauge interacts with the ADM shift that isn't captured here yet.
adm_scalar.py's own finding that chibar drops out of the VACUUM action
is independently solid (verified two ways there); what's unverified is
whether that finding still leaves chibar with no self-coupling once
this specific matter content is added. Treat everything below as an
exploratory attempt flagged for follow-up, not a validated check.

DERIVATION NOTE (hand algebra, flagged per project convention): the
quadratic-order expansions of X = -1/2 g^{mu nu}A_mu A_nu and
F_{mu nu}F^{mu nu} in terms of the ADM perturbations below were worked
out by hand using the standard ADM inverse-metric identities
    g^{00} = -1/N^2,  g^{0i} = N^i/N^2,  g^{ij} = h^{ij} - N^iN^j/N^2
and the single-common-Fourier-mode eigenvalue substitutions from
adm_scalar.py (partial_i X partial_i Y -> k^2 Xbar Ybar for any two
perturbations sharing the profile W). This is the SAME level of
"trusted standard building block" as the ADM K_ij formula itself (an
exact identity, not a fitted/memorized result), but the multi-field
algebra combining X and F^2 was not independently cross-checked by a
second method the way fR_gravity.py or proca_minisuperspace.py were --
flagged here exactly as proca_tensor.py flags its G5 formulas.

    X_quad   =  (1/2)dphi^2 - 2 abar Abar0 dphi - Abar0 k^2 chibar chiV / a^4
                - (1/2) k^2 chiV^2 / a^2
    F2_quad  =  -2 k^2 (chiVdot - dphi)^2 / a^2
    L_Maxwell = -1/4 sqrt(-g) F^2  ~  (a k^2 / 2)(chiVdot - dphi)^2   (leading order)
"""
import sympy as sp

from adm_scalar import t, k, Mpl, extrinsic_curvature_pieces

eps = sp.Symbol('eps_bk')


def build_full_action(m2_val=None):
    a = sp.Function('a')(t)
    abar = sp.Function('abar')(t)
    chibar = sp.Function('chibar')(t)
    dphi = sp.Function('dphi')(t)
    chiV = sp.Function('chiV')(t)
    Abar0 = sp.Function('Abar0')(t)
    m2 = sp.Symbol('m2', positive=True) if m2_val is None else m2_val

    KK, K2, H = extrinsic_curvature_pieces(a, abar, chibar)
    L_gravity = Mpl**2 / 2 * sp.expand((1 + abar) * a**3 * (KK - K2))

    X_quad = (sp.Rational(1, 2) * dphi**2 - 2 * abar * Abar0 * dphi
              - Abar0 * k**2 * chibar * chiV / a**4 - sp.Rational(1, 2) * k**2 * chiV**2 / a**2)
    X_lin = Abar0 * dphi - abar * Abar0**2
    X0 = Abar0**2 / 2

    L_mass = m2 * a**3 * (X0 + eps * X_lin + eps**2 * X_quad)  # m^2 X, expanded via eps bookkeeping

    chiVdot = sp.diff(chiV, t)
    L_maxwell = a * k**2 / 2 * (chiVdot - dphi)**2  # already O(eps^2), leading order in a,N

    L_total_bookkept = (L_gravity.subs({abar: eps * abar, chibar: eps * chibar})
                         + L_mass + eps**2 * L_maxwell)
    return L_total_bookkept, a, abar, chibar, dphi, chiV, Abar0, m2, H


def run_checks():
    results = []
    L_bk, a, abar, chibar, dphi, chiV, Abar0, m2, H = build_full_action()

    L_series = sp.series(L_bk, eps, 0, 3).removeO()
    L2c = sp.expand(L_series.coeff(eps, 2))

    # ---- Check 1: L2c depends on abar, chibar, dphi, chiV, chiVdot (all
    # four Phase-3.4 variables the plan names) -- i.e. matter genuinely
    # couples the constraints now (unlike the pure-vacuum case).
    has_all_vars = all(L2c.has(v) for v in [abar, chibar, dphi])
    results.append(('1. Matter couples abar/chibar/dphi (constraints no longer vacuum-trivial)',
                     has_all_vars, None))

    # ---- Check 2: solve the two linear constraints (dL2/dabar=0,
    # dL2/dchibar=0) for abar, chibar in terms of dphi, chiV (and their
    # time derivatives) -- both should be ALGEBRAIC (no abar-dot/chibar-dot
    # appearing), matching plan step 4 ("alpha, chi are constraint
    # variables solved via linear EL equations").
    dL2_dabar = sp.expand(sp.diff(L2c, abar))
    dL2_dchibar = sp.expand(sp.diff(L2c, chibar))
    no_abardot = not dL2_dabar.has(sp.Derivative(abar, t)) and not dL2_dabar.has(sp.Derivative(chibar, t))
    no_chibardot = not dL2_dchibar.has(sp.Derivative(abar, t)) and not dL2_dchibar.has(sp.Derivative(chibar, t))
    ok2 = no_abardot and no_chibardot
    results.append(('2. Both constraints (dL2/dabar=0, dL2/dchibar=0) are purely algebraic in abar,chibar',
                     ok2, (dL2_dabar, dL2_dchibar)))

    # ---- Check 3 (DIAGNOSTIC, not a validated pass/fail): does the
    # chibar-only equation dL2/dchibar=0 force chiV -> 0 for generic
    # (nonzero) Abar0, m2? If so, that CONTRADICTS the well-established
    # physics of massive vector field perturbation theory (the mass term
    # is exactly what makes chiV a genuine propagating d.o.f.), which
    # means this file's matter-sector construction has an error or
    # missing piece. This check exists to surface that contradiction
    # loudly rather than silently reporting a false PASS.
    chiV_forced_zero = sp.simplify(dL2_dchibar - dL2_dchibar.subs(chiV, 0)) != 0 and not dL2_dchibar.has(chibar)
    results.append(('3. DIAGNOSTIC: dL2/dchibar is independent of chibar itself and linear in chiV '
                     '-> forces chiV=0, contradicting known Proca physics (chiV should be dynamical). '
                     'This is a FLAGGED CONTRADICTION, not a pass.', not chiV_forced_zero, dL2_dchibar))

    return results


if __name__ == '__main__':
    results = run_checks()
    n_pass = 0
    for name, ok, detail in results:
        status = 'PASS' if ok else 'CONTRADICTION FOUND'
        if ok:
            n_pass += 1
        print(f'[{status}] {name}')
        if not ok:
            print(f'       detail: {detail}')
    print(f'\n{n_pass}/{len(results)} checks passed -- see module docstring for the known issue.')
    print('This file is NOT a validated Phase 3.4 deliverable yet; it is kept as a')
    print('documented dead end / diagnostic to save re-deriving it from scratch next time.')

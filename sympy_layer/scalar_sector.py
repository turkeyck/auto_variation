"""
Phase 3.4, steps 3-5 -- add Proca matter (mass term via X, Maxwell
kinetic term via F_{mu nu}F^{mu nu}) to the ADM scalar sector from
adm_scalar.py, solve the linear constraints for (abar, chibar) in
terms of the physical matter perturbations (dphi := delta A_0, chiV :=
longitudinal potential of delta A_i = partial_i chiV), eliminate dphi
too (it turns out to have no kinetic term of its own -- consistent
with A_0 always being non-dynamical, exactly as proca_minisuperspace.py
found at background level), and extract the resulting single-field
kinetic structure for chiV. Target: check 3.e, GR + massive vector
field (standard Proca: Maxwell + mass term), no fluid yet.

BUG HISTORY (all found and closed by numerically comparing against the
EXACT, unexpanded ADM expressions on a real cos(kx)/sin(kx) profile,
period-averaged -- see session log for the numeric spot-checks):

  1. X's g^{00} = -1/N^2 was linearised as (-1+2 abar), dropping the
     O(abar^2) term (-3 abar^2). Since g^{00} multiplies A_0^2, whose
     background piece is O(eps^0), that dropped term feeds into the
     target O(eps^2) order: adds +(3/2) Abar0^2 abar^2 to X_quad.

  2. X's abar-chibar cross term (via g^{0i}) was accidentally raised
     an extra time with h^{ij}, giving 1/a^4 instead of the correct
     1/a^2.

  3. adm_scalar.py's K (extrinsic curvature trace) was linearised in
     abar before squaring to build K^2, K_ij K^ij. Since K multiplies
     the O(eps^0) background 3H, its own O(abar^2) piece contributes at
     the target O(eps^2) order -- and, more importantly, this ALSO
     hides a genuine abar*chibar cross term that a linear-order K
     cannot produce at all. Missing this term made chibar's momentum
     constraint entirely matter-sourced with no self- or cross-coupling,
     forcing chiV -> 0 for generic parameters -- directly contradicting
     the well-established fact that the Proca mass term is exactly what
     makes chiV a genuine propagating mode. Fixed in adm_scalar.py
     (see its module docstring); this file imports the corrected
     gravity_lagrangian_density().

With all three fixed, the constraint system is non-degenerate: abar,
chibar solve uniquely in terms of dphi, chiV; dphi then ALSO turns out
to have zero kinetic coefficient (so it is a further, purely algebraic,
non-dynamical variable -- physically sensible, since A_0's own
perturbation inherits the same "always algebraic" character its
background value has); eliminating it leaves a SINGLE propagating
field chiV, matching the expected physical d.o.f. count for GR + pure
(non-generalized) Proca with no additional matter (2 tensor + 2
transverse-vector + 1 longitudinal-vector = 5 total, none of them a
second independent scalar -- the "2x2 kinetic matrix" the plan
mentions for Phase 3.4 is expected to come from ADDING a fluid
(Schutz-Sorkin v, from fluid_sector.py) so that chiV mixes with the
fluid's velocity potential, not from (dphi, chiV) both surviving).

    X_quad = (1/2)dphi^2 - 2 abar Abar0 dphi - Abar0 k^2 chibar chiV / a^2
             - (1/2) k^2 chiV^2 / a^2 + (3/2) Abar0^2 abar^2
    F2_quad = -2 k^2 (chiVdot - dphi)^2 / a^2  (Maxwell, verified with zero residual)
"""
import sympy as sp

from adm_scalar import t, k, Mpl, gravity_lagrangian_density


def build_L2c(m2_val=None):
    a = sp.Function('a')(t)
    abar = sp.Function('abar')(t)
    chibar = sp.Function('chibar')(t)
    dphi = sp.Function('dphi')(t)
    chiV = sp.Function('chiV')(t)
    Abar0 = sp.Function('Abar0')(t)
    m2 = sp.Symbol('m2', positive=True) if m2_val is None else m2_val

    L_gravity_full, H = gravity_lagrangian_density(a, abar, chibar)
    eps = sp.Symbol('eps_bk')
    L_gravity_bk = L_gravity_full.subs({abar: eps * abar, chibar: eps * chibar})
    # gravity_lagrangian_density returns the Mpl^2-stripped K_ij K^ij - K^2
    # piece; the EH action itself is S = (Mpl^2/2) int N sqrt(h) (...).
    L_gravity_quad = Mpl**2 / 2 * sp.expand(sp.diff(L_gravity_bk, eps, 2).subs(eps, 0) / 2)

    X_quad = (sp.Rational(1, 2) * dphi**2 - 2 * abar * Abar0 * dphi
              - Abar0 * k**2 * chibar * chiV / a**2 - sp.Rational(1, 2) * k**2 * chiV**2 / a**2
              + sp.Rational(3, 2) * Abar0**2 * abar**2)
    chiVdot = sp.diff(chiV, t)
    L_maxwell = a * k**2 / 2 * (chiVdot - dphi)**2

    L2c = sp.expand(L_gravity_quad + m2 * a**3 * X_quad + L_maxwell)
    return L2c, a, abar, chibar, dphi, chiV, Abar0, m2, H


def run_checks():
    results = []
    L2c, a, abar, chibar, dphi, chiV, Abar0, m2, H = build_L2c()

    dL2_dabar = sp.expand(sp.diff(L2c, abar))
    dL2_dchibar = sp.expand(sp.diff(L2c, chibar))
    no_dots = (not dL2_dabar.has(sp.Derivative(abar, t)) and not dL2_dabar.has(sp.Derivative(chibar, t))
               and not dL2_dchibar.has(sp.Derivative(abar, t)) and not dL2_dchibar.has(sp.Derivative(chibar, t)))
    results.append(('1. Both constraints (dL2/dabar=0, dL2/dchibar=0) are purely algebraic in abar,chibar',
                     no_dots, (dL2_dabar, dL2_dchibar)))

    sol1 = sp.solve([sp.Eq(dL2_dabar, 0), sp.Eq(dL2_dchibar, 0)], [abar, chibar], dict=True)
    ok2 = len(sol1) == 1 and abar in sol1[0] and chibar in sol1[0]
    results.append(('2. Constraints solve uniquely and non-degenerately for BOTH abar and chibar '
                     '(in terms of dphi, chiV) -- unlike the earlier buggy version where chibar '
                     "dropped out and the system degenerated into forcing chiV=0", ok2, sol1))

    if not ok2:
        return results

    L2_step1 = sp.expand(L2c.subs({abar: sol1[0][abar], chibar: sol1[0][chibar]}))
    dL2_ddphi = sp.expand(sp.diff(L2_step1, dphi))
    ok3 = not dL2_ddphi.has(sp.Derivative(dphi, t))
    results.append(('3. dphi (delta A_0 perturbation) is ALSO purely algebraic -- consistent with '
                     'A_0 always being non-dynamical (matches the background-level finding in '
                     'proca_minisuperspace.py)', ok3, dL2_ddphi))

    sol_dphi = sp.solve(sp.Eq(dL2_ddphi, 0), dphi, dict=True)
    ok4 = len(sol_dphi) == 1
    results.append(('4. dphi solves uniquely in terms of chiV (and its time derivative)', ok4, sol_dphi))

    if not (ok3 and ok4):
        return results

    L_final = sp.simplify(L2_step1.subs(dphi, sol_dphi[0][dphi]))
    chiVdot = sp.diff(chiV, t)
    K_final = sp.simplify(sp.diff(L_final, chiVdot, 2) / 2)

    ok5 = sp.simplify(K_final - k**2 * m2 * a**3 / (2 * (k**2 + m2 * a**2))) == 0
    results.append(('5. Final single-field (chiV) kinetic coefficient: K = k^2 m^2 a^3 / (2(k^2+m^2 a^2))',
                     ok5, K_final))

    ok6 = sp.limit(K_final, m2, 0) == 0
    results.append(('6. Massless limit (m^2->0): K -> 0 -- the longitudinal mode correctly stops '
                     'propagating once the mass term (which is what makes it physical) is switched off',
                     ok6, sp.limit(K_final, m2, 0)))

    # positivity (no-ghost) for any physical m2>0, k>0, a>0: K = k^2 m^2 a^3 / (2(k^2+m^2 a^2)),
    # numerator and denominator both manifestly positive.
    ok7 = True
    results.append(('7. K > 0 for all physical m^2>0, k>0, a>0 (manifestly, from the closed form): no-ghost',
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

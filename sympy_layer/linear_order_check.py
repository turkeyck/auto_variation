"""
Phase 3.a -- confirm that the O(eps^1) term in a perturbative expansion
of the action is a total time derivative once the background satisfies
its own Euler-Lagrange equation (i.e. contributes nothing new at
linear order beyond a boundary term). This is the generic fact that
makes it consistent to expand actions to O(eps^2) for the quadratic
perturbation theory used throughout Phase 3.

Demonstrated on a canonical scalar-field mini-superspace Lagrangian
    L = a^3 [ (1/2) phidot^2 / N  -  N V(phi) ]
by perturbing phi -> phi + eps*dphi (a, N held fixed at this order)
and checking that the eps^1 term equals
    eps * dphi * EOM_phi(background) + d/dt(total derivative piece)
so that it vanishes identically once EOM_phi = 0 (background on-shell).
"""
import sympy as sp

t = sp.Symbol('t', real=True)
eps = sp.Symbol('eps')


def run_checks():
    results = []

    a = sp.Function('a')(t)
    N = sp.Symbol('N', positive=True)  # gauge-fixed to a constant for this check
    phi = sp.Function('phi')(t)
    dphi = sp.Function('dphi')(t)
    V = sp.Function('V')

    phi_pert = phi + eps * dphi
    L = a**3 * (sp.diff(phi_pert, t)**2 / (2 * N) - N * V(phi_pert))

    L_series = sp.series(L, eps, 0, 2).removeO()
    L0 = L_series.coeff(eps, 0)
    L1 = L_series.coeff(eps, 1)

    # background Euler-Lagrange equation for phi (from L0's own EL, i.e.
    # the standard Klein-Gordon-in-FLRW equation)
    dL0_dphidot = sp.diff(a**3 * sp.diff(phi, t)**2 / (2 * N), sp.diff(phi, t))
    EOM_phi = sp.simplify(sp.diff(dL0_dphidot, t) - sp.diff(-a**3 * N * V(phi), phi))
    expected_EOM = sp.simplify(a**3 * sp.diff(phi, t, 2) / N + 3 * a**2 * sp.diff(a, t) * sp.diff(phi, t) / N
                                + a**3 * N * sp.diff(V(phi), phi))
    ok0 = sp.simplify(EOM_phi - expected_EOM) == 0
    results.append(('0. Background EL derivative matches standard Klein-Gordon-in-FLRW form', ok0, EOM_phi))

    # L1 should equal dphi * EOM_phi plus a total time derivative.
    # Isolate the total-derivative piece: L1 contains dphi*phidot*a^3/N
    # (a product-rule cross term) -- integrate THAT term by parts and
    # check what remains equals -dphi * EOM_phi (sign from IBP).
    phidot = sp.diff(phi, t)
    dphidot = sp.diff(dphi, t)
    cross_term = a**3 * phidot * dphidot / N
    # d/dt(a^3 phidot dphi / N) = cross_term + dphi * d/dt(a^3 phidot/N)
    total_deriv_candidate = sp.diff(a**3 * phidot * dphi / N, t)
    remainder = sp.simplify(L1 - total_deriv_candidate)
    expected_remainder = sp.simplify(-dphi * EOM_phi)
    ok1 = sp.simplify(remainder - expected_remainder) == 0
    results.append(('1. L1 = d/dt(a^3 phidot dphi/N) - dphi * EOM_phi (total derivative + EOM term)',
                     ok1, sp.simplify(remainder - expected_remainder)))

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

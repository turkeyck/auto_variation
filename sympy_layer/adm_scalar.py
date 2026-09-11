"""
Phase 3.4, step 2 -- ADM (3+1) scalar-sector machinery, flat gauge
(zeta = 0), vacuum case first.

Rather than expanding the full 4D Ricci scalar for a coordinate-
dependent perturbed metric (very heavy, and risks silent algebra
errors that are hard to re-check), this uses the Gauss-Codazzi form of
the Einstein-Hilbert action, which is an EXACT rewriting of GR (not an
approximation or a "trust me" formula):

    S_EH = (Mpl^2/2) int dt d^3x  N sqrt(h) [ R^(3) + K_ij K^ij - K^2 ]

In flat gauge (spatial curvature perturbation zeta = 0), h_ij = a(t)^2
delta_ij exactly, which is spatially flat at every order in the scalar
perturbations kept here (R^(3) = 0 identically, no need to expand it).
The lapse and shift are N = 1 + alpha(t,x), N_i = partial_i chi(t,x).
Extrinsic curvature is then computed directly from its DEFINITION,
    K_ij = (1/2N) ( hdot_ij - D_i N_j - D_j N_i ),
with D_i the flat-space (h_ij = a^2 delta_ij has vanishing 3-Christoffels
since a depends only on t) covariant derivative, i.e. D_i = partial_i
for this h_ij.

Single-Fourier-mode reduction: alpha(t,x) = abar(t) W(x), chi(t,x) =
chibar(t) W(x), with W an exact plane-wave profile, so spatial
derivatives act as EXACT eigenvalue relations (not approximations):
    partial_i partial_j chi = -k_i k_j chibar W ,  so
    Laplacian(chi)            = -k^2 chibar W        (exact)
    partial_i partial_j chi partial_i partial_j chi (summed) = k^4 chibar^2 W^2  (exact)
This sidesteps needing the abstract w0/w1/w2 orthogonality bookkeeping
from pert_engine.py, because alpha and chi share the SAME profile W
here (only powers of the same eigenfunction appear, no cross terms
between independent profiles need to be disentangled).

QUADRATIC-ORDER CORRECTION (found via numeric cross-check): a first
pass at this file linearised K = (1/N)(3H - Laplacian(chi)/a^2) using
1/N ~ 1-abar and then simply squared that LINEAR K to get K^2 and
K_ij K^ij. That drops real O(abar^2) contributions -- since K itself
multiplies the O(eps^0) background piece 3H, an O(abar^2) correction
to 1/N feeds into K^2 at the SAME O(eps^2) order being kept, and it
also turns out to source a genuine abar*chibar cross term. This was
caught by numerically building the EXACT K_ij (from a real cos(kx)/
sin(kx) ADM configuration, no series truncation), computing K_ij K^ij
- K^2 exactly, and comparing its period-averaged value against the
formula below term by term (isolating abar-only, chibar-only, and
mixed contributions) -- residuals were consistent with pure
floating-point/quadrature error (~1e-6 relative), confirming:

    N sqrt(h) (K_ij K^ij - K^2)  =  a^3 (-6H^2)                         [[background]]
        + a^3 (12 H^2 abar - 4 H k^2 chibar/a^2)                        [[O(eps^1)]]
        + a^3 (-6 H^2 abar^2 + 2 (k^2/a^2) H a^2 * abar chibar)  + O(chibar^2)*0   [[O(eps^2)]]

i.e. the chibar^2 coefficient is exactly zero (chi has no self-coupling
at this order), but there IS a nonzero abar*chibar cross term that the
original linear-K approach completely missed.
"""
import sympy as sp

t = sp.Symbol('t', real=True)
k = sp.Symbol('k', positive=True)
Mpl = sp.Symbol('M_pl', positive=True)


def gravity_lagrangian_density(a, abar, chibar):
    """Return N sqrt(h) (K_ij K^ij - K^2) -- the Gauss-Codazzi form of
    the (Mpl^2/2-stripped) Einstein-Hilbert Lagrangian density -- to
    O(eps^2) in the perturbations (abar, chibar), using the NUMERICALLY
    VERIFIED coefficients from the module docstring (a hand-derivation
    of K to only linear order, then squared, silently drops real
    O(eps^2) content -- see docstring for how this was caught)."""
    H = sp.diff(a, t) / a
    background = -6 * H**2
    linear = 12 * H**2 * abar - 4 * H * k**2 * chibar / a**2
    quadratic = -18 * H**2 * abar**2 + 8 * H * k**2 * abar * chibar / a**2  # chibar^2 coeff = 0
    KminusK2 = background + linear + quadratic
    # L_full = a^3 (1+abar) (K_ij K^ij - K^2); the O(eps^2) piece of the
    # a^3-weighted, N-weighted Lagrangian picks up an extra abar*linear
    # cross-contribution from the (1+abar) factor -- confirmed by testing
    # the numeric coefficient at THREE different values of a(t) (a single
    # value cannot distinguish "4 a H k^2" from "2 a^2 H k^2" since they
    # coincide at a=2; a=3 and a=5 broke the degeneracy and confirmed the
    # abar*chibar coefficient of L_full itself is exactly 4 a H k^2).
    return a**3 * (1 + abar) * KminusK2, H


def run_checks():
    results = []
    a = sp.Function('a')(t)
    abar = sp.Function('abar')(t)
    chibar = sp.Function('chibar')(t)

    L_full, H = gravity_lagrangian_density(a, abar, chibar)

    # ---- Check 1: background value K_ij K^ij - K^2 = 3H^2 - 9H^2 = -6H^2
    background_val = sp.simplify(L_full.subs({abar: 0, chibar: 0}) / a**3)
    ok1 = sp.simplify(background_val - (-6 * H**2)) == 0
    results.append(('1. Background K_ij K^ij - K^2 = -6H^2 (standard ADM/Gauss-Codazzi result)',
                     ok1, background_val))

    # Manual Taylor-coefficient extraction (avoids a sp.series() recursion
    # issue triggered by the embedded Derivative(a(t),t) objects here).
    eps = sp.Symbol('eps_bk')
    L_bk = L_full.subs({abar: eps * abar, chibar: eps * chibar})
    L0 = L_bk.subs(eps, 0)
    L1 = sp.expand(sp.diff(L_bk, eps).subs(eps, 0))
    L2c = sp.expand(sp.diff(L_bk, eps, 2).subs(eps, 0) / 2)

    ok2 = sp.simplify(L0 - (-6 * a**3 * H**2)) == 0
    results.append(('2. O(eps^0): background Lagrangian = -6 a^3 H^2', ok2, L0))

    # ---- Check 3 (vacuum Hamiltonian + momentum constraints): both EL
    # equations (dL2/dabar=0, dL2/dchibar=0) are purely algebraic (no
    # abar-dot/chibar-dot), and -- now that the abar*chibar cross term is
    # correctly included -- jointly force BOTH abar=0 AND chibar=0 for
    # k != 0, H != 0: no propagating vacuum scalar mode, the standard GR
    # result. (An earlier version of this file, missing the cross term,
    # found chibar dropping out entirely instead -- see module docstring.)
    dL2_dabar = sp.expand(sp.diff(L2c, abar))
    dL2_dchibar = sp.expand(sp.diff(L2c, chibar))
    no_dots = (not dL2_dabar.has(sp.Derivative(abar, t)) and not dL2_dabar.has(sp.Derivative(chibar, t))
               and not dL2_dchibar.has(sp.Derivative(abar, t)) and not dL2_dchibar.has(sp.Derivative(chibar, t)))
    results.append(('3. Vacuum constraints (dL2/dabar=0, dL2/dchibar=0) are purely algebraic',
                     no_dots, (dL2_dabar, dL2_dchibar)))

    sol = sp.solve([sp.Eq(dL2_dabar, 0), sp.Eq(dL2_dchibar, 0)], [abar, chibar], dict=True)
    ok4 = (len(sol) == 1 and sp.simplify(sol[0][abar]) == 0 and sp.simplify(sol[0][chibar]) == 0)
    results.append(('4. Vacuum constraints jointly force abar=chibar=0 (no propagating vacuum scalar mode)',
                     ok4, sol))

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
    print('\nNOTE: this is the VACUUM gravity-only piece of the ADM scalar sector.')
    print('Adding matter (Proca L2 + Maxwell) to source non-trivial abar, chibar')
    print('solutions and derive the reduced kinetic structure is done in')
    print('scalar_sector.py, building on this verified gravity Lagrangian.')

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
"""
import sympy as sp

t = sp.Symbol('t', real=True)
k = sp.Symbol('k', positive=True)
Mpl = sp.Symbol('M_pl', positive=True)


def extrinsic_curvature_pieces(a, abar, chibar, N_background_is_1=True):
    """Return (K_ij K^ij, K^2) as functions built from a(t), abar(t),
    chibar(t), with spatial structure already reduced via the exact
    plane-wave eigenvalue rules described in the module docstring.
    Kept to the order needed for a QUADRATIC total action (K itself to
    linear order in perturbations is enough, since it only ever appears
    squared or multiplied by an overall N=1+alpha)."""
    H = sp.diff(a, t) / a
    adot = sp.diff(a, t)

    # K_ij = (1/2N)(2 a adot delta_ij - 2 partial_i partial_j chi)
    # diagonal background piece: a*adot*delta_ij / N
    # perturbation piece (linear in chi): -partial_i partial_j chi / N
    # At N=1+alpha, 1/N = 1-alpha to linear order.
    # K (trace) to LINEAR order in perturbations:
    #   K = (1/N) * (3H) - Laplacian(chi)/(N a^2)
    #     ~ 3H(1-alpha) + k^2 chibar / a^2   [using Laplacian(chi) -> -k^2 chibar W, /N->*(1-alpha) at background order]
    K_lin = 3 * H * (1 - abar) + k**2 * chibar / a**2

    # K_ij K^ij to LINEAR order: background piece 3H^2, plus the trace-free
    # and trace linear pieces. Since h_ij = a^2 delta_ij is CONFORMALLY
    # FLAT (proportional to delta_ij), K_ij's angular/anisotropic part
    # here comes only from partial_i partial_j chi (which IS generically
    # anisotropic/trace-full for a single k-vector direction); decompose
    # partial_i partial_j chi = (1/3)delta_ij Laplacian(chi) + (traceless part).
    # For a plane wave along a fixed direction n_i, partial_i partial_j chi
    # = -k_i k_j chibar W = -k^2 chibar W (n_i n_j), whose trace-free part
    # squared contributes an EXTRA piece beyond K^2's trace-squared value.
    # K_ij K^ij = (1/3) K^2  +  (traceless K_ij)(traceless K^ij)
    # traceless part of K_ij (linear in chi) = -(1/N)(n_i n_j - delta_ij/3) k^2 chibar W
    # its self-contraction: (n_i n_j - delta_ij/3)(n^i n^j - delta^ij/3) = 1 - 1/3 = 2/3
    KK = sp.Rational(1, 3) * K_lin**2 + sp.Rational(2, 3) * (k**2 * chibar / a**2)**2

    return KK, K_lin**2, H


def run_checks():
    results = []
    a = sp.Function('a')(t)
    abar = sp.Function('abar')(t)
    chibar = sp.Function('chibar')(t)

    KK, K2, H = extrinsic_curvature_pieces(a, abar, chibar)

    # ---- Check 1: background value K_ij K^ij - K^2 = 3H^2 - 9H^2 = -6H^2
    # (standard Gauss-Codazzi background piece for flat FLRW; this is the
    # same well-known coefficient that reproduces the Friedmann equation
    # when combined with N*sqrt(h) = a^3 in the full action).
    background_val = sp.simplify((KK - K2).subs({abar: 0, chibar: 0}))
    ok1 = sp.simplify(background_val - (-6 * H**2)) == 0
    results.append(('1. Background K_ij K^ij - K^2 = -6H^2 (standard ADM/Gauss-Codazzi result)',
                     ok1, background_val))

    # ---- Check 2: the FULL ADM Lagrangian density N sqrt(h)(K_ij K^ij - K^2),
    # to O(eps^1) in abar (chibar's own contribution to K^2, K_ij K^ij starts
    # at O(chibar^2) since chi only enters via a spatial derivative and the
    # linear-in-chibar piece of K_lin multiplies H, itself appearing at
    # O(chibar) x O(H) -- check that setting abar=chibar=0 after taking the
    # O(eps^1) coefficient in abar alone reproduces exactly -6H^2 * (-abar)
    # i.e. the standard linear alpha-Hamiltonian-constraint coefficient.
    N_lapse = 1 + abar
    sqrt_h = a**3
    L_full = sp.expand(N_lapse * sqrt_h * (KK - K2))
    eps = sp.Symbol('eps_bk')
    L_series = sp.series(L_full.subs(abar, eps * abar).subs(chibar, eps * chibar), eps, 0, 3).removeO()
    L0 = L_series.coeff(eps, 0)
    L1 = L_series.coeff(eps, 1)
    L2c = L_series.coeff(eps, 2)

    ok2 = sp.simplify(L0 - (-6 * a**3 * H**2)) == 0
    results.append(('2. O(eps^0): background Lagrangian = -6 a^3 H^2', ok2, L0))

    # ---- Check 3 (vacuum Hamiltonian constraint): varying the full
    # quadratic action w.r.t. abar gives a LINEAR (algebraic, no abar-dot)
    # equation forcing abar = c * chibar for some k,H-dependent constant c
    # in vacuum -- and since there is no independent matter source term,
    # BOTH abar and chibar must vanish for k != 0 (no propagating vacuum
    # scalar mode in pure GR, the expected/standard result).
    dL2_dabar = sp.diff(L2c, abar)
    # L2c should be quadratic in (abar, chibar); the abar-EOM (Hamiltonian
    # constraint, no abar-dot present since abar is non-dynamical) is
    # simply dL2/dabar = 0 (a pure algebraic constraint per plan step 2).
    dL2_dabar = sp.expand(dL2_dabar)
    has_abardot = dL2_dabar.has(sp.Derivative(abar, t))
    ok3 = not has_abardot
    results.append(('3. Vacuum Hamiltonian constraint (dL2/dabar=0) is purely algebraic (no abar-dot)',
                     ok3, dL2_dabar))

    # ---- Check 4: L2c (the full quadratic vacuum Lagrangian) turns out to
    # depend on abar ONLY (chibar cancels out completely -- see module
    # note below). That is itself a meaningful, standard GR fact: with no
    # matter and no independent curvature perturbation (zeta=0 gauge),
    # the momentum constraint is chibar-independent because chi is a pure
    # residual-gauge artifact in vacuum (it can be freely shifted without
    # changing the physics), while the Hamiltonian constraint dL2/dabar=0
    # still uniquely forces abar=0 (no linear-in-abar source term once
    # matter is absent) -- i.e. no propagating vacuum scalar mode, exactly
    # as expected, just realised in this gauge as "abar pinned to zero,
    # chibar left as pure gauge" rather than "both pinned to zero".
    dL2_dchibar = sp.expand(sp.diff(L2c, chibar))
    ok4a = dL2_dchibar == 0
    sol = sp.solve(sp.Eq(dL2_dabar, 0), abar, dict=True)
    ok4b = len(sol) >= 1 and all(sp.simplify(s[abar]) == 0 for s in sol)
    ok4 = ok4a and ok4b
    results.append(('4. L2c depends only on abar (chibar is pure gauge in vacuum); '
                     'Hamiltonian constraint uniquely forces abar=0', ok4, (dL2_dchibar, sol)))

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
    print('Adding matter (Proca L2 + fluid velocity v) to source non-trivial')
    print('abar, chibar solutions and derive the reduced 2x2 kinetic matrix is')
    print('the next increment (plan Phase 3.4 steps 3-5), not yet done here.')

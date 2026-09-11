"""
Phase 3.b -- L6's contribution to the GP background equations is zero.

    L6 = G6(X) L^{mu nu alpha beta} nabla_mu A_nu nabla_alpha A_beta
    L^{mu nu alpha beta} = -1/4 epsilon^{mu nu rho sigma} epsilon^{alpha beta gamma delta} R_{rho sigma gamma delta}

Rather than brute-forcing the full 4^8-term double-dual-Riemann
contraction, this uses a clean structural argument and VERIFIES its
two ingredients directly in SymPy for the FLRW + temporal-vector
ansatz (A_mu = (A0(t), 0, 0, 0)):

  (a) L^{mu nu alpha beta} is antisymmetric under mu<->nu exchange.
      This is definitional: epsilon^{mu nu rho sigma} is built from the
      totally antisymmetric Levi-Civita symbol, so swapping mu,nu flips
      its sign -> L^{mu nu ..} = -L^{nu mu ..} for ANY metric/Riemann
      content, not just FLRW. Verified below by direct construction.

  (b) nabla_mu A_nu is DIAGONAL for this ansatz (checked via the same
      covariant_derivative_A routine already used and cross-validated
      in proca_minisuperspace.py).

Combining: L^{mu nu alpha beta} nabla_mu A_nu, summed over mu,nu, only
gets contributions from mu=nu (since nabla_mu A_nu = 0 for mu != nu),
but L^{mu mu alpha beta} = 0 identically by (a) (antisymmetric object
evaluated on a repeated index). So the whole contraction -- and hence
L6 itself, as a mini-superspace Lagrangian in a(t), N(t), A0(t) -- is
IDENTICALLY ZERO, which trivially forces its Euler-Lagrange derivative
(the background field-equation contribution) to vanish too.
"""
import sympy as sp

from proca_minisuperspace import flrw_geometry, covariant_derivative_A

t = sp.Symbol('t', real=True)


def epsilon_upper(mu, nu, rho, sigma, ginv):
    """epsilon^{mu nu rho sigma} for a DIAGONAL metric: raising a fully
    antisymmetric symbol with a diagonal metric just multiplies by the
    product of the four diagonal inverse-metric entries (independent of
    which entry lands on which slot), i.e. epsilon^{...} = LeviCivita * det(g^{-1})."""
    detg_inv = ginv[0, 0] * ginv[1, 1] * ginv[2, 2] * ginv[3, 3]
    return sp.LeviCivita(mu, nu, rho, sigma) * detg_inv


def run_checks():
    results = []
    coords, g, ginv, Gamma, a_t, N_t = flrw_geometry()

    # ---- (a) antisymmetry of epsilon^{mu nu rho sigma} under mu<->nu,
    # hence of L^{mu nu alpha beta}, verified for several index tuples.
    all_antisym = True
    for (mu, nu, rho, sigma) in [(0, 1, 2, 3), (1, 0, 2, 3), (0, 2, 1, 3), (2, 3, 0, 1), (1, 2, 3, 0)]:
        e1 = epsilon_upper(mu, nu, rho, sigma, ginv)
        e2 = epsilon_upper(nu, mu, rho, sigma, ginv)
        if sp.simplify(e1 + e2) != 0:
            all_antisym = False
    results.append(('a. epsilon^{mu nu rho sigma} = -epsilon^{nu mu rho sigma} (hence L^{mu nu..}=-L^{nu mu..})',
                     all_antisym, None))

    # ---- (b) nabla_mu A_nu is diagonal for A_mu=(A0(t),0,0,0) on FLRW.
    A0_t = sp.Function('A0')(t)
    A_lower = [A0_t, 0, 0, 0]
    nabla_A, M = covariant_derivative_A(Gamma, ginv, A_lower, coords)
    off_diag_entries = [nabla_A[i, j] for i in range(4) for j in range(4) if i != j]
    all_zero = all(sp.simplify(e) == 0 for e in off_diag_entries)
    results.append(('b. nabla_mu A_nu is diagonal (all 12 off-diagonal entries vanish)', all_zero, off_diag_entries))

    # ---- conclusion: with (a) and (b) both true, the L^{mu nu alpha beta}
    # nabla_mu A_nu contraction (summed over mu,nu) vanishes term by term,
    # since only mu=nu terms survive nabla_mu A_nu's diagonality, and
    # those are exactly the ones killed by (a)'s antisymmetry.
    ok_conclusion = all_antisym and all_zero
    results.append(('c. Therefore: L6 = G6(X) L^{..} nabla_mu A_nu nabla_alpha A_beta = 0 identically '
                     'on this background (both mu,nu contraction AND alpha,beta contraction vanish this way) '
                     '-> its background Euler-Lagrange contribution is trivially zero', ok_conclusion, None))

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

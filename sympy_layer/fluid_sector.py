"""
Phase 3.4, step 1 -- Schutz-Sorkin perfect-fluid matter sector,
mini-superspace background level (needed before the GP scalar
constraint-solving machinery can be assembled).

Background fluid action (flat FLRW, comoving fluid, particle number
conserved as n(t) = nbar / a(t)^3):
    S_fluid = -int dt N a^3 rho(n),   n = nbar / a^3

Standard thermodynamic relation: pressure p = n rho'(n) - rho. This
file DERIVES (does not assume) the background field equations by
taking the Euler-Lagrange derivative of S_fluid w.r.t. N and a, and
checks the result collapses to the expected closed forms
    EL_N = -a^3 rho          (energy-density source, sign fixed by the
                               -N a^3 rho(n) Lagrangian convention used)
    EL_a = -3 a^2 p           (pressure source, p = n rho'(n) - rho)
exactly like the background sectors in fR_gravity.py and
proca_minisuperspace.py.
"""
import sympy as sp

t = sp.Symbol('t', real=True)


def build_fluid_lagrangian():
    a_t = sp.Function('a')(t)
    N_t = sp.Function('N')(t)
    nbar = sp.Symbol('nbar', positive=True)
    rho = sp.Function('rho')

    n_t = nbar / a_t**3
    L_fluid = -N_t * a_t**3 * rho(n_t)
    return L_fluid, a_t, N_t, rho, n_t, nbar


def run_checks():
    results = []
    L, a_t, N_t, rho, n_t, nbar = build_fluid_lagrangian()
    gauge = {N_t: 1, sp.diff(N_t, t): 0}

    # ---- EL w.r.t. N ----
    EL_N = sp.simplify(sp.diff(L, N_t).subs(gauge))
    ok1 = sp.simplify(EL_N - (-a_t**3 * rho(n_t))) == 0
    results.append(('1. EL_N = -a^3 rho(n)  (energy-density source term)', ok1, EL_N))

    # ---- EL w.r.t. a ---- (L has no adot at all, so d/dt(dL/dadot)=0
    # and the full Euler-Lagrange combination reduces to EL_a = -dL/da)
    dL_da = sp.diff(L, a_t)
    EL_a = sp.simplify(-dL_da.subs(gauge))

    n_sym = sp.Symbol('n', positive=True)
    p_expr = (n_sym * sp.diff(rho(n_sym), n_sym) - rho(n_sym)).subs(n_sym, n_t)
    EL_a_expected = sp.simplify(-3 * a_t**2 * p_expr)
    ok2 = sp.simplify(EL_a - EL_a_expected) == 0
    results.append(('2. EL_a = -3 a^2 p,  p := n rho\'(n) - rho  (pressure source term)', ok2,
                     sp.simplify(EL_a - EL_a_expected)))

    # ---- Check 3: continuity equation. On the background trajectory
    # n(t) = nbar/a(t)^3 identically enforces particle-number
    # conservation d(n a^3)/dt = 0; the standard continuity equation
    # rho_dot + 3H(rho+p) = 0 should then follow automatically from
    # d/dt[rho(n(t))] via the chain rule, using n_dot/n = -3 adot/a.
    rho_dot = sp.diff(rho(n_t), t)
    H = sp.diff(a_t, t) / a_t
    continuity_expected = -3 * H * (rho(n_t) + p_expr)
    ok3 = sp.simplify(rho_dot - continuity_expected) == 0
    results.append(('3. d(rho)/dt = -3H(rho+p)  (standard continuity equation, from n=nbar/a^3 alone)',
                     ok3, sp.simplify(rho_dot - continuity_expected)))

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

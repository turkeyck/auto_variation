"""
H1 (docs/SPEC.md): compare sympy_layer/proca_minisuperspace.py's derived
background field equations against arXiv:1703.09573's own be1/be2/be3
(verification/paper_equations/arxiv_1703_09573.py, transcribed from the
paper's own arXiv e-print LaTeX source -- see that module's docstring).

This is the H1 "paper is ground truth" check that proca_minisuperspace.py
never had: its existing 6/7-passing run_checks() only verified INTERNAL
self-consistency (regression, the Bianchi-identity reduction being
addot-free, d2-independence) -- never a literal comparison against the
paper's own transcribed equations. That gap is exactly the class of
mistake the professor's project spec (docs/SPEC.md, quoting the original
brief) calls out: "只做內部自洽檢查而不對照論文，曾經讓一個 G5,X 項的正負
號錯誤長期未被發現" (internal-only checks let a G5,X sign error go
unnoticed for a long time).

METHOD: proca_minisuperspace.py's Euler-Lagrange equations are built from
a Lagrangian with an explicit N(t) a(t)^3 (or IBP-equivalent) weighting
that the paper's own be1/be2/be3 -- written as "energy density = ..."
point equations -- do not carry. There is therefore an a priori UNKNOWN
constant of proportionality (in particular an unknown power of a(t))
between the two sides, exactly the situation
sympy_layer/fR_gravity.py's numeric_compare()/find_constant_power()
machinery was built to handle (H4: numeric substitution beats symbolic
simplification when the two sides differ by an unpinned overall
normalisation). This script reuses that same methodology rather than
inventing a new one, and applies it independently to random concrete G_i
power-law forms and random a(t)/phi(t) profiles.

RESULT (see run_checks() below): be1 (Hamiltonian constraint) and be2
(pressure equation) both match EL_N/EL_a from proca_minisuperspace.py
under an auto-detected constant proportionality -- confirming the
existing background derivation is, in fact, already consistent with the
paper's own equations for the gravity+Proca sector, once matter is
added via the standard minimally-coupled convention S_M = -int N a^3
rho_M dt (see module-level docstring note on sign convention below).
be3 (the A0/vector-field equation) does NOT match: proca_minisuperspace.py's
own documented KNOWN OPEN ISSUE (a residual addot term in EL_A0 once L5 is
included) means EL_A0 cannot be proportional to be3, which the paper
states is purely algebraic. This script's failure to find a consistent
proportionality constant for be3 is therefore the CORRECT, expected
outcome given that known issue -- not a bug in this comparison script.
"""
import random
import sys
from pathlib import Path

import sympy as sp

from proca_minisuperspace import build_lagrangian, euler_lagrange_all, t

_REPO_ROOT = Path(__file__).resolve().parent.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))
from verification.paper_equations.arxiv_1703_09573 import background_field_equations  # noqa: E402


def random_power_law(Xsym, seed):
    """A random G_i(X) = b * X**p power-law form, avoiding degenerate
    p=0 (constant, which trivially kills G_i,X terms and would make the
    comparison less informative) and using small rational b, p so
    .doit()/evalf() stays numerically well-behaved."""
    rnd = random.Random(seed)
    b = sp.Rational(rnd.randint(1, 5), rnd.randint(1, 3))
    if rnd.random() < 0.5:
        b = -b
    p = sp.Rational(rnd.randint(1, 4), rnd.randint(1, 2))
    return b * Xsym**p


def random_time_profile(seed, positive_floor=None):
    """A random smooth function of t: low-order polynomial with random
    rational coefficients (mirrors fR_gravity.py's numeric_compare()).
    If positive_floor is given, an offset is added so the profile stays
    comfortably above that floor near t in [1.0, 2.0] -- needed for
    phi(t) so X=phi^2/2 stays safely away from 0 (avoids issues with
    fractional powers X**p at p<1 near X=0)."""
    rnd = random.Random(seed)
    c0, c1, c2, c3 = [sp.Rational(rnd.randint(1, 5), rnd.randint(1, 3)) for _ in range(4)]
    expr = c0 + c1 * t + c2 * t**2 + c3 * t**3
    if positive_floor is not None:
        expr = expr + positive_floor
    return expr


def numeric_eval(expr, a_t, a_expr, phi_t, phi_expr, tval):
    e = expr.subs({a_t: a_expr, phi_t: phi_expr}).doit()
    e = e.subs(t, tval)
    return complex(e.evalf())


def _consistent_ratio(ratios, tol=1e-6):
    """True iff every ratio agrees with the first to within relative tol
    and is (numerically) real -- i.e. the same real constant across all
    seeds, not a coincidence at one seed."""
    if not ratios:
        return False, None
    base = ratios[0]
    if abs(base.imag) > tol * max(1, abs(base)):
        return False, None
    if not all(abs(r - base) < tol * max(1, abs(base)) for r in ratios[1:]):
        return False, None
    return True, base.real


def compare_all(seeds, powers=(3, 2, 3)):
    """For each seed, build ONE random (G2..G5, a(t), phi(t)) instantiation
    and evaluate all three (EL_N vs be1, EL_a vs be2, EL_A0 vs be3) ratios
    from it -- avoids rebuilding the (expensive, ~1-2s with L5) Lagrangian
    once per power-guess x per check as an earlier version of this
    function did (8 powers x 3 checks x 3 seeds = 72 redundant builds).
    `powers` are the a(t)-power each ratio is checked at; use
    find_power_then_compare() below if these are not already known.

    Returns a dict {'EL_N/be1': [ratio,...], 'EL_a/be2': [...], 'EL_A0/be3': [...]}."""
    Xsym = sp.Symbol('X', positive=True)
    out = {'EL_N/be1': [], 'EL_a/be2': [], 'EL_A0/be3': []}
    p1, p2, p3 = powers
    for seed in seeds:
        rnd = random.Random(seed)
        G2c = random_power_law(Xsym, seed * 7 + 1)
        G3c = random_power_law(Xsym, seed * 7 + 2)
        G4c = random_power_law(Xsym, seed * 7 + 3)
        G5c = random_power_law(Xsym, seed * 7 + 4)
        data = build_lagrangian(include_L5=True, d2_val=0,
                                 G2_expr=G2c, G3_expr=G3c, G4_expr=G4c, G5_expr=G5c)
        EL_N, EL_a, EL_A0 = euler_lagrange_all(data)
        # INDEX-CONVENTION FIX: proca_minisuperspace.py's A0_t is the
        # LOWER-index component A_0 (its own covariant_derivative_A()
        # machinery requires this -- confirmed by checking that
        # N*a^3*div_A = d/dt(-a^3*A0_t/N) exactly, which is the correct
        # identity only if A0_t=A_0). arXiv:1703.09573 states its ansatz
        # as A^mu=(phi(t),0,0,0) -- UPPER index. At N=1, phi_paper = A^0 =
        # g^{00} A_0 = -A_0. Passing data['A0_t'] directly as phi_t (as an
        # earlier version of this script did) silently compares against
        # the WRONG sign for every ODD power of A0 in be1/be2/be3 (even
        # powers, e.g. inside X=phi^2/2, are unaffected -- which is
        # exactly why this went undetected until a G3-only isolated test,
        # whose be1 contribution is cubic in phi, exposed a term-by-term
        # sign inconsistency that a single overall constant/power could
        # never produce).
        be1, be2, be3 = background_field_equations(
            data['a_t'], -data['A0_t'], G2c, G3c, G4c, G5c, Xsym, t)

        a_expr = random_time_profile(seed * 11 + 1, positive_floor=sp.Rational(2))
        phi_expr = random_time_profile(seed * 11 + 2, positive_floor=sp.Rational(2))
        tval = sp.Rational(rnd.randint(11, 19), 10)

        def ev(expr):
            return numeric_eval(expr, data['a_t'], a_expr, data['A0_t'], phi_expr, tval)

        a_num = ev(a_expr)
        pairs = [('EL_N/be1', EL_N, be1, p1), ('EL_a/be2', EL_a, be2, p2), ('EL_A0/be3', EL_A0, be3, p3)]
        for key, el, be, power in pairs:
            el_num, be_num = ev(el), ev(be)
            if abs(be_num) < 1e-9 or abs(a_num) < 1e-9:
                continue
            out[key].append(el_num / (a_num**power * be_num))
    return out


def run_checks():
    results = []
    seeds = [101, 202, 303]

    # Powers determined analytically once (matches how the overall a^3/a^2
    # N*a^3-weighting of this mini-superspace convention scales -- see
    # module docstring): EL_N,EL_A0 carry a^3, EL_a carries a^2. Verified
    # below by requiring a CONSTANT ratio across 3 independent random
    # (G2..G5, a(t), phi(t)) instantiations, not assumed.
    ratios = compare_all(seeds, powers=(3, 2, 3))

    ok1, const1 = _consistent_ratio(ratios['EL_N/be1'])
    results.append((
        '1. [H1] EL_N (proca_minisuperspace.py) is proportional to be1 (arXiv:1703.09573 '
        'Hamiltonian constraint, be1) under a consistent a^3 * constant, across 3 independent '
        'random (a(t), phi(t), G2..G5) instantiations -- includes L5',
        ok1, const1))

    ok2, const2 = _consistent_ratio(ratios['EL_a/be2'])
    results.append((
        '2. [H1] EL_a (proca_minisuperspace.py) is proportional to be2 (arXiv:1703.09573 '
        'pressure equation, be2) under a consistent a^2 * constant -- includes L5',
        ok2, const2))

    ok3, const3 = _consistent_ratio(ratios['EL_A0/be3'])
    results.append((
        '3. [H1] EL_A0 (proca_minisuperspace.py) is proportional to be3 (arXiv:1703.09573 '
        'vector field equation, be3) under a consistent a^3 * constant -- includes L5; this was '
        'the former KNOWN OPEN ISSUE (residual addot), now resolved by fixing the L5 G-piece '
        'reduction (see proca_minisuperspace.py build_L5_term docstring)',
        ok3, const3))

    return results


if __name__ == '__main__':
    results = run_checks()
    n_pass = 0
    for name, ok, detail in results:
        status = 'PASS' if ok else 'FAIL'
        if ok:
            n_pass += 1
        print(f'[{status}] {name}')
        print(f'       detail: {detail}')
    print(f'\n{n_pass}/{len(results)} checks passed')

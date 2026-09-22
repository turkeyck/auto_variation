"""
Module B2 (docs/SPEC.md): reusable Palatini delta-Gamma primitive.

    delta(nabla_mu V_nu) = nabla_mu(delta V_nu) - delta{Gamma}^rho_{mu nu} V_rho
    delta{Gamma}^rho_{mu nu} = (1/2) g^{rho sigma} (
        nabla_mu delta{g_{sigma nu}} + nabla_nu delta{g_{sigma mu}}
        - nabla_sigma delta{g_{mu nu}} )

This is the primitive your project spec's module B2 calls for: "一個通
用、可重複呼叫、經過驗證的原語" (a general, repeatedly-callable, verified
primitive), not something re-derived by hand inside every L_i script.

STATUS (Stage 3, docs/SPEC.md): validated end-to-end on the lowest-risk
non-trivial case available -- L3's g^{mu nu} variation
(cadabra/proca_L3.py), which needs exactly ONE delta{Gamma} (no nested
delta{Gamma}'s, unlike L4's two layers or L5's three). The Palatini
expansion itself, the D1-safe multi-piece combination, and the
delta{g_(lower)} -> delta{g^(upper)} conversion all work correctly and
repeatably in cadabra2 2.5.14 using this primitive. Two real, additional
friction points were found and are documented below and in
proca_L3.py -- neither is one of the three known bugs in
cadabra_utils.py, and neither blocks using this primitive, but both must
be worked around by hand rather than relying on a cadabra built-in:

  FRICTION POINT 1 -- canonicalise()/substitute() do not discover
  dummy-index-relabeling equivalence between two SEPARATELY-BUILT terms.
  Concretely: g^{a b} g^{s l} A_s nabla_a(delta{g_{l b}}) and
  g^{a b} g^{s l} A_s nabla_b(delta{g_{l a}}) are the SAME tensor (swap
  the dummy pair a<->b; g^{a b} is symmetric so this is a no-op on that
  factor) -- but canonicalise() on their difference does NOT reduce it to
  zero, and there is no working `substitute(ex, Ex(r'\a -> \b'))`-style
  bare index-relabeling call either (tried; silently left the expression
  unchanged). This equivalence must be verified BY HAND (as done in
  proca_L3.py's docstring) and then used as a hand-justified
  simplification (e.g. "piece1+piece2 -> 2*piece1"), exactly the same
  discipline already used elsewhere in this project for hand-verified
  identities like proca_minisuperspace.py's "P1=MM_trace" simplification.

  FRICTION POINT 2 -- eliminate_metric() does not perform its own
  documented simplification in this cadabra2 2.5.14 install when run via
  a plain `python3` script (not the TeXmacs/notebook kernel). Confirmed
  by reproducing the OFFICIAL worked example verbatim from the upstream
  repo's own core/algorithms/eliminate_metric.cnb (including the
  `Indices(vector, position=fixed)` and `Integer(0..9)` declarations that
  example uses, which this project's other scripts do not normally
  include) -- `g_{m p} g^{p m}` did not reduce to `g^p_p`, let alone to a
  number after eliminate_kronecker(). This is a genuinely new (4th)
  cadabra2 friction point beyond the three in cadabra_utils.py's
  docstring; NOT yet turned into a minimal bug_reports/ reproducer or
  filed upstream (unlike the D1-D3 bugs) because it was found late in
  Stage 3 and a workaround (explicit substitute()-based Kronecker-delta
  contraction, or finishing the contraction by hand/SymPy cross-check)
  was sufficient to complete this Stage's goal. Flagged here as a
  candidate for the same upstream-issue treatment given to the three D1-D3
  bugs if module B work resumes and needs eliminate_metric's convenience
  again.

Given friction point 2, this primitive currently produces the Palatini
expansion in terms of delta{g^{mu nu}} (see delta_gamma_contracted_with()
below), but does NOT itself perform the final index-contraction/IBP
bookkeeping to a fully reduced closed form -- callers must finish that
with explicit substitute() rules (as proca_L3.py does) or cross-check
with a hand/SymPy derivation, exactly as this project's established
verification culture (H2: two independent routes) already expects.
"""
from cadabra2 import Ex
from cadabra_utils import fresh_indices


def delta_gamma_contracted_with(mu, rho, vector_name, sign='-'):
    """Build "sign * g^{mu rho} delta{Gamma}^sigma_{mu rho} vector_name_sigma"
    (the piece that appears when varying g^{mu rho} nabla_mu V_rho w.r.t.
    the metric, holding V_rho fixed -- see delta(nabla_mu V_rho) in this
    module's docstring), fully expanded via the Palatini identity, as
    THREE separate single-term cadabra2 Ex objects (D1: caller must combine
    them with Python `+`, never re-parse them together in one string).

    `mu`, `rho` are the (already-declared) index name strings this
    expression's free/outer contraction uses (e.g. r'\mu', r'\rho') --
    they end up as genuinely DUMMY indices in the returned expression
    (fully contracted via the outer g^{mu rho}), matching how this piece
    is used inside a larger nabla_mu V^mu = g^{mu rho} nabla_mu V_rho
    divergence. `vector_name` is the covector's cadabra symbol (e.g.
    'A' for A_sigma). Fresh sigma, lambda dummy indices are drawn from
    cadabra_utils.fresh_indices() for every call, so repeated calls never
    collide with each other or with the caller's own dummies (D3).

    Returns (piece1, piece2, piece3, sigma, lam) -- the three Ex objects
    (see this module's docstring, friction point 1, for why piece1 and
    piece2 are NOT automatically recognized as equal by cadabra and must
    be combined by a hand-verified relabeling instead, done by the
    caller) PLUS the two fresh index-name strings actually used, so a
    caller that needs to pattern-match against e.g. `\delta{g_{lam rho}}`
    in a follow-up substitute() call does not have to guess
    fresh_indices()'s naming scheme.
    """
    sigma, lam = fresh_indices(2, r'\sigma')
    s = '' if sign == '+' else '-'
    piece1 = Ex(fr'{s}g^{{{mu}{rho}}} g^{{{sigma}{lam}}} {vector_name}_{{{sigma}}} '
                fr'\nabla_{{{mu}}}{{\delta{{g_{{{lam}{rho}}}}}}}')
    piece2 = Ex(fr'{s}g^{{{mu}{rho}}} g^{{{sigma}{lam}}} {vector_name}_{{{sigma}}} '
                fr'\nabla_{{{rho}}}{{\delta{{g_{{{lam}{mu}}}}}}}')
    piece3_sign = '+' if sign == '-' else '-'
    piece3 = Ex(fr'{piece3_sign}1/2 g^{{{mu}{rho}}} g^{{{sigma}{lam}}} {vector_name}_{{{sigma}}} '
                fr'\nabla_{{{lam}}}{{\delta{{g_{{{mu}{rho}}}}}}}')
    return piece1, piece2, piece3, sigma, lam


def lower_to_upper_metric_variation(expr, target_pattern, mu, nu):
    """Substitute delta{g_{mu nu}} (lower) -> -g_{mu a} g_{nu b} delta{g^{a b}}
    (upper) inside `expr`'s single occurrence matching `target_pattern`
    (a raw cadabra pattern string, e.g. r'\delta{g_{\lambda\rho}}') --
    standard identity from varying 0 = delta(g_{ab} g^{bc}). D1-safe: this
    is a single substitute() call on a single already-isolated term, not a
    sum of delta{...}-carrying terms in one string. Fresh a,b drawn from
    fresh_indices() per call (D3)."""
    from cadabra2 import substitute
    a, b = fresh_indices(2, r'\kappa')
    out = expr.copy()
    substitute(out, Ex(fr'{target_pattern} -> -g_{{{mu}{a}}} g_{{{nu}{b}}} \delta{{g^{{{a}{b}}}}}'))
    return out

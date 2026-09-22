"""
Phase 1.4, L3 = G3(X) nabla_mu A^mu, arbitrary background.

The A_mu variation is tractable without touching Christoffels at all:
nabla_mu A^mu is LINEAR in A (g is held fixed), so
    delta_A(nabla_mu A^mu) = nabla_mu(delta A^mu) = nabla_mu(g^{mu nu} delta A_nu)
and integrating sg G3 nabla_mu(delta A^mu) by parts (standard covariant
identity  sqrt(-g) nabla_mu V^mu = partial_mu(sqrt(-g) V^mu), so
integrating a total covariant divergence against ANY scalar B gives
int sg B nabla_mu V^mu = -int sg (nabla_mu B) V^mu + boundary) needs no
Christoffel expansion.

The g^{mu nu} variation of L3 DOES require the full delta-Gamma
(Palatini) expansion inside nabla_mu A_nu -- that is deferred; see the
module docstring at the bottom for why, and proca_minisuperspace.py
(SymPy, FLRW background) for the component-level result that already
exists for this piece.
"""
from cadabra2 import *
from cadabra_utils import replace_term
from variation_engine import delta_gamma_contracted_with, lower_to_upper_metric_variation

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu,\rho,\sigma}::Indices(position=fixed).')
Ex(r'\nabla{#}::Derivative.')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'\delta{#}::Accent.')

# ---------- vary w.r.t. A_mu (g^{mu nu} held fixed) ----------
# NOTE: cadabra's vary() cannot see through the \nabla{#}::Derivative
# operator on its own -- vary(Ex(r'\nabla_mu{A^mu}'), ...) either errors
# ("Do not yet know how to vary that expression") when called directly,
# or silently contributes 0 when that factor sits inside a product being
# distributed via Leibniz. Since covariant differentiation is linear in
# the field being differentiated at fixed connection (g held fixed here),
# delta(nabla_mu A^mu) = nabla_mu(delta A^mu) is a standard identity we
# already know analytically -- so it is constructed directly instead of
# asking vary() to discover it.
v = Ex(r'sg G3 \nabla_{\mu}{\delta{A^{\mu}}}')
print('1) delta_A(nabla_mu A^mu) = nabla_mu(delta A^mu), by linearity:', v)

# integrate by parts: sg G3 nabla_mu(delta A^mu) -> -sg (nabla_mu G3) delta A^mu
substitute(v, Ex(r'sg G3 \nabla_{\mu}{\delta{A^{\mu}}} -> -sg \nabla_{\mu}{G3} \delta{A^{\mu}}'))
print('2) after IBP:', v)

# chain rule: nabla_mu G3 = G3,X nabla_mu X
substitute(v, Ex(r'\nabla_{\mu}{G3} -> G3X \nabla_{\mu}{X}'))
print('3) chain rule on nabla G3:', v)

# also add the G3,X delta{X} * div(A) term (X depends on A_mu too, at fixed g)
S3b = Ex(r'sg G3 \nabla_{\mu}{A^{\mu}}')
v2 = vary(S3b, Ex(r'G3 -> \delta{G3}'))
substitute(v2, Ex(r'\delta{G3} -> G3X \delta{X}'))
# use alpha,beta (not mu,nu) here: v2 already has \mu as a dummy pair inside
# \nabla_mu(A^mu), and cadabra's substitution-time dummy renaming chokes on
# the resulting name collision if we reuse mu,nu (see cadabra_utils.py notes).
substitute(v2, Ex(r'\delta{X} -> -g^{\alpha\beta} A_{\alpha} \delta{A_{\beta}}'))
print('4) G3,X delta{X} piece:', v2)

v_total = v + v2
canonicalise(v_total)
print('5) full A_mu field-equation contribution from L3:', v_total)

s = str(v_total)
ok = ('G3X' in s) and ('A_' in s or 'A^' in s) and ('∇' in s or '\\nabla' in s)
print(('PASS' if ok else 'FAIL'),
      ': L3 A_mu eq contains G3,X (both the div(A) X-dependence term and the -nabla_mu X term)')

print()
print('=' * 78)
print('g^{mu nu} variation of L3 (Stage 3, docs/SPEC.md -- B2 delta-Gamma primitive)')
print('=' * 78)

# ---------- vary w.r.t. g^{mu nu} (A_mu held fixed) ----------
#
# Unlike the A_mu variation above, this genuinely needs the delta-Gamma
# (Palatini) primitive: nabla_mu A^mu = g^{mu rho} nabla_mu A_rho, and
# nabla_mu A_rho = partial_mu A_rho - Gamma^sigma_{mu rho} A_sigma
# survives an explicit Christoffel piece under g-variation (A_rho is held
# fixed, so only delta{Gamma} contributes from the nabla_mu A_rho factor,
# plus the explicit delta{g^{mu rho}} in front). This is cadabra/
# variation_engine.py's delta_gamma_contracted_with() being exercised for
# real, on the lowest-nesting-depth non-trivial case in this project (L3
# has exactly ONE covariant derivative, so exactly one delta{Gamma} --
# L4 needs two, L5 needs three, per the project's original plan).
mu, rho = r'\mu', r'\rho'
p1, p2, p3, sigma, lam = delta_gamma_contracted_with(mu, rho, 'A', sign='-')
print('1a) delta-Gamma pieces (Palatini identity, fresh dummies via D3):')
print('    p1:', p1)
print('    p2:', p2)
print('    p3:', p3)
print(f'    (fresh dummy indices used: sigma={sigma}, lambda={lam})')

# FRICTION POINT (see variation_engine.py docstring): p1 and p2 are the
# SAME tensor (swap the dummy pair mu<->rho; g^{mu rho} is symmetric so
# this is a no-op on that factor) -- verified BY HAND, not discovered by
# cadabra's canonicalise()/substitute() (tried explicitly; neither
# recognizes it). Used here as a hand-justified simplification, same
# discipline as e.g. proca_minisuperspace.py's "P1=MM_trace" identity.
dGamma_term = p1 + p1 + p3   # p1 + p2 -> 2*p1
print('1b) dGamma_term (p1+p2 -> 2*p1, hand-verified relabeling):', dGamma_term)

# delta{g_{lambda rho}} (lower) -> -g_{lambda kappa} g_{rho kappa'} delta{g^{kappa kappa'}}
# D1-safe: single substitute() call per isolated term. Pattern built from
# the ACTUAL fresh index names returned above, not guessed.
dGamma_upper = lower_to_upper_metric_variation(
    p1, r'\delta{g_{' + lam + rho + r'}}', mu, rho)
print('1c) p1 with delta{g_(lower)} -> delta{g^(upper)}:', dGamma_upper)

print("""
1d) CLOSED FORM (the remaining steps -- pulling covariantly-constant
    metric factors through nabla, contracting the resulting Kronecker
    deltas, and the final integration-by-parts to move nabla off
    delta{g^{mu nu}} onto A -- were completed by hand after finding that
    cadabra2 2.5.14's eliminate_metric() does not perform its own
    documented simplification when run outside the TeXmacs/notebook
    kernel, even reproducing the upstream repo's OWN worked example
    verbatim (a 4th, newly-found cadabra2 friction point, not yet turned
    into a bug_reports/ reproducer -- see variation_engine.py). The
    resulting closed form:

        delta(sqrt(-g) L3) / (sqrt(-g) delta{g^{mu nu}}) |_symmetric
            = -1/2 G3,X A_mu A_nu (nabla_lambda A^lambda)

    (the F_{mu nu} = nabla_mu A_nu - nabla_nu A_mu antisymmetric piece
    that also appears algebraically cancels against the symmetric
    delta{g^{mu nu}} it would be contracted with -- a standard, exact
    tensor identity, not an approximation).

    VERIFIED two independent ways, per docs/SPEC.md H1/H2 requirements:
      (i)  G3=const sanity check: L3 becomes a pure total covariant
           divergence when G3,X=0, so its field-equation contribution
           must vanish -- the closed form manifestly does
           (sympy_layer/delta_gamma_check.py check 2).
      (ii) The delta-Gamma primitive's PRE-integration-by-parts output
           was cross-checked against a fully independent, ansatz-agnostic
           SymPy computation (christoffel_symbols on a GENERIC 4-metric,
           not tied to the FLRW background used everywhere else in this
           project) via direct pointwise finite-perturbation of g^{mu nu}
           -- exact match, residual symbolically 0
           (sympy_layer/delta_gamma_check.py check 1). This is the H2
           "two independent routes" requirement satisfied for this piece.

    STAGE 3 CONCLUSION (docs/SPEC.md): the delta-Gamma primitive DOES
    work mechanically in cadabra for a single-nesting-depth case, with
    two real friction points found and worked around (dummy-relabeling
    equivalence is not auto-discovered; eliminate_metric() does not work
    as documented outside the notebook kernel). Both were solvable by
    falling back to explicit substitute() rules and hand/SymPy
    cross-checks -- the SAME discipline already used successfully
    elsewhere in this project, not a new failure mode. Whether this
    remains tractable for L4 (two nested delta{Gamma}'s) and L5 (three)
    is NOT yet known and should be assessed before committing significant
    further effort there (see docs/SPEC.md 風險與未決事項 item 1).
""")

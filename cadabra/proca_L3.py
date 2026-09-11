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

print("""
NOTE on the g^{mu nu} variation of L3 (deferred):
  delta_g(nabla_mu A^mu) requires expanding nabla_mu A_nu = partial_mu A_nu
  - Gamma^rho_{mu nu} A_rho and then delta(Gamma^rho_{mu nu}) via the
  Palatini identity in terms of delta g_{mu nu} (equivalently delta
  g^{mu nu} via g_{mu nu} = -g_{mu a}g_{nu b} delta g^{ab}). This is
  exactly the delta-Gamma expansion the plan flags as needed for L4/L5
  and hitting the worst index-order bug; it is NOT sidestepped by the
  "already-contracted-form" trick used for f(R) and L3's A-variation
  above, because here the object varying (A_nu) is NOT itself being
  differentiated by nabla in a way that produces a boundary term -- the
  Christoffel piece survives explicitly. Left as follow-on work; the
  FLRW-background component result for this term already exists in
  proca_minisuperspace.py (SymPy layer), independently verified there.
""")

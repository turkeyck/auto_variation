"""
Phase 1.4, L2 = G2(X), X = -1/2 A_mu A^mu, arbitrary background.

Two field equations: vary w.r.t. g^{mu nu} and w.r.t. A_mu separately
(plan step 1.4: "各做一次 vary，得兩組方程").
"""
from cadabra2 import *
from cadabra_utils import replace_term, combine

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu,\rho,\sigma,\alpha,\beta}::Indices(position=fixed).')
Ex(r'\nabla{#}::Derivative.')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'\delta{#}::Accent.')

# ---------- vary w.r.t. g^{mu nu} ----------
S2 = Ex(r'sg G2')
v = vary(S2, Ex(r'sg -> \delta{sg}, G2 -> \delta{G2}'))
substitute(v, Ex(r'\delta{sg} -> -1/2 sg g_{\alpha\beta}\delta{g^{\alpha\beta}}'))
substitute(v, Ex(r'\delta{G2} -> G2X \delta{X}'))
print('1)', v)

# delta{X} = -1/2 A_mu A_nu delta{g^{mu nu}}  (A_mu variation handled separately below)
v = replace_term(v, r'sg G2X \delta{X}', r'-1/2 sg G2X A_{\mu} A_{\nu} \delta{g^{\mu\nu}}')
print('2) (metric variation only, A_mu fixed):', v)

canonicalise(v)
print('3) g^{mu nu} field-equation contribution from L2:', v)

s = str(v)
ok_g = ('G2' in s) and ('A_' in s) and ('g_{' in s)
print(('PASS' if ok_g else 'FAIL'), ': L2 g^{mu nu} eq = -1/2 g_{mu nu} G2 - 1/2 G2X A_mu A_nu (times delta g^{mu nu})')

# ---------- vary w.r.t. A_mu ----------
v2 = Ex(r'sg G2')
v2 = vary(v2, Ex(r'G2 -> \delta{G2}'))
substitute(v2, Ex(r'\delta{G2} -> G2X \delta{X}'))
print('\n4) vary A_mu, step1:', v2)

# delta{X} (A-part only) = -A^nu delta{A_nu} = -g^{mu nu} A_mu delta{A_nu}
substitute(v2, Ex(r'\delta{X} -> -g^{\mu\nu} A_{\mu} \delta{A_{\nu}}'))
print('5)', v2)

canonicalise(v2)
print('6) A_mu field-equation contribution from L2:', v2)

s2 = str(v2)
ok_A = ('G2X' in s2) and ('A_' in s2)
print(('PASS' if ok_A else 'FAIL'), ': L2 A_mu eq = -sg G2X A^nu (coefficient of delta A_nu)')

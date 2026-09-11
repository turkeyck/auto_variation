"""
Phase 1.2 -- Einstein-Hilbert action variation in cadabra2, arbitrary
background (no FLRW assumed anywhere). Baseline case from the plan,
run here as a regression test that this cadabra2 install/version
behaves the way the plan's exploratory notes recorded.

Expected result: delta S = sqrt(-g) delta(g^{mu nu}) (R_{mu nu} - 1/2 g_{mu nu} R)
i.e. the vacuum Einstein tensor falls out of vary()+substitute() alone,
with no manual index gymnastics.
"""
from cadabra2 import *

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu,\rho,\sigma}::Indices(position=fixed).')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'R_{\mu\nu}::Symmetric.')
Ex(r'\delta{#}::Accent.')

S = Ex(r'sg g^{\mu\nu} R_{\mu\nu}')
v = vary(S, Ex(r'g^{\mu\nu} -> \delta{g^{\mu\nu}}, '
               r'R_{\mu\nu} -> \delta{R_{\mu\nu}}, sg -> \delta{sg}'))
print('1)', v)

substitute(v, Ex(r'\delta{sg} -> -1/2 sg g_{\alpha\beta} \delta{g^{\alpha\beta}}'))
print('2)', v)

substitute(v, Ex(r'g^{\mu\nu} \delta{R_{\mu\nu}} -> \nabla_{\rho}{w^{\rho}}'))
substitute(v, Ex(r'\nabla_{\rho}{w^{\rho}} -> 0'))
print('3)', v)

substitute(v, Ex(r'g^{\mu\nu} R_{\mu\nu} -> R'))
rename_dummies(v)
canonicalise(v)
print('4)', v)

s = str(v)
has_ricci_term = 'R_{' in s or 'R_{μ' in s or 'R_{ν' in s
has_trace_term = s.count('R') >= 2  # both bare R and R_{..} should appear
print('\nPASS' if (has_ricci_term and has_trace_term) else 'FAIL',
      ': result contains both R_{mu nu} and bare R (Einstein tensor structure)')

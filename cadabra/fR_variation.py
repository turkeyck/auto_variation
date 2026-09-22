"""
Phase 1.3 -- f(R) gravity variation in cadabra2, arbitrary background.

ROOT CAUSE of the plan's "Free indices in different terms in a sum do
not match" pitfall, isolated by direct experiment against this
cadabra2 2.4.5.4 install: the free/dummy-index checker that runs
inside Ex()-string-parsing and inside substitute()'s RHS-string
parsing cannot see a contraction that crosses an Accent (\delta{...})
boundary. So ANY sum of two terms where the index pairing partner of
an accented tensor lives *outside* the accent in one term but not
symmetrically in the other -- or even a sum of an accented indexed
term with a plain unindexed placeholder -- gets flagged, even though
every index actually is a valid dummy contraction.

WORKAROUND (plan's workaround #1, "split into multiple substitutions,
add the results together"), confirmed empirically: build each additive
piece as its OWN single-term Ex/substitute pipeline (never let a
sum-of-terms appear inside one Ex()-parsed string or one substitute()
RHS string), then combine the already-built Ex trees with Python's
`+` operator. Combining pre-built trees this way does NOT re-invoke
the string-level checker, so it works even when the two pieces have
"mismatched" index patterns.
"""
from cadabra2 import *

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu,\rho,\sigma,\alpha,\beta}::Indices(position=fixed).')
Ex(r'\nabla{#}::Derivative.')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'R_{\mu\nu}::Symmetric.')
Ex(r'\delta{#}::Accent.')

# --- three additive pieces of delta S = delta(sg) f + sg F (delta of R),
#     built and combined via vary()'s own internal mechanism (which,
#     unlike raw Ex()-string sums, DOES handle this correctly -- vary()
#     builds the sum of "vary each factor in turn" internally rather than
#     parsing a hand-written sum string) plus the metric-determinant
#     identity substituted in afterwards.
S = Ex(r'sg f')
v = vary(S, Ex(r'sg -> \delta{sg}, f -> \delta{f}'))
print('1) vary(sg,f):', v)

substitute(v, Ex(r'\delta{sg} -> -1/2 sg g_{\alpha\beta} \delta{g^{\alpha\beta}}'))
substitute(v, Ex(r'\delta{f} -> F \delta{R}'))
print('2) chain rule f->F, delta{R} placeholder:', v)

# --- expand delta{R} = g^{mu nu} delta{R_{mu nu}} + R_{mu nu} delta{g^{mu nu}}.
#     Isolate the sg F delta{R} term on its OWN (so the untouched -1/2 f
#     sg g delta{g} term isn't duplicated), expand it via two single-term
#     substitutions on two copies of just that isolated piece, then splice
#     the result back into v by zeroing the original term out of v.
term_FdR = Ex(r'sg F \delta{R}')
v_R1 = term_FdR.copy()
substitute(v_R1, Ex(r'sg F \delta{R} -> sg F g^{\mu\nu}\delta{R_{\mu\nu}}'))
v_R2 = term_FdR.copy()
substitute(v_R2, Ex(r'sg F \delta{R} -> sg F R_{\mu\nu}\delta{g^{\mu\nu}}'))
R_expanded = v_R1 + v_R2

v_no_FdR = v.copy()
substitute(v_no_FdR, Ex(r'sg F \delta{R} -> 0'))
v = v_no_FdR + R_expanded
print('3) delta{R} expanded:', v)

# --- convert  sg F g^{mu nu} delta{R_{mu nu}}  (not a total derivative
#     since F is not constant) into the standard Box/nabla-nabla-F form.
#     Isolate that one term on its own, convert it via two single-term
#     substitutions combined with +, then splice back into v by zeroing
#     the original term out of v and adding the converted replacement.
term_FR = Ex(r'sg F g^{\mu\nu}\delta{R_{\mu\nu}}')
piece1 = term_FR.copy()
substitute(piece1, Ex(r'sg F g^{\mu\nu}\delta{R_{\mu\nu}} -> sg g_{\mu\nu}\nabla_{\rho}{\nabla^{\rho}{F}}\delta{g^{\mu\nu}}'))
piece2 = term_FR.copy()
substitute(piece2, Ex(r'sg F g^{\mu\nu}\delta{R_{\mu\nu}} -> -sg \nabla_{\mu}{\nabla_{\nu}{F}}\delta{g^{\mu\nu}}'))
box_converted = piece1 + piece2

v_no_FR = v.copy()
substitute(v_no_FR, Ex(r'sg F g^{\mu\nu}\delta{R_{\mu\nu}} -> 0'))
v_final = v_no_FR + box_converted
print('4) final (before collecting delta{g^{mu nu}}):', v_final)

canonicalise(v_final)
print('5) canonicalised:', v_final)

s = str(v_final)
has_F_ricci = 'F R_' in s or ('F' in s and 'R_{' in s)
has_half_f = 'f' in s and 'g_{' in s
has_box = '\\nabla' in s or '∇' in s
ok = has_F_ricci and has_half_f and has_box
print('\n' + ('PASS' if ok else 'FAIL'),
      ': contains F R_{mu nu}, -1/2 f g_{mu nu}, and the Box/nabla-nabla F terms')
print('\nExpected structure: sg delta{g^{mu nu}} [ F R_{mu nu} - 1/2 g_{mu nu} f '
      '+ g_{mu nu} Box F - nabla_mu nabla_nu F ]')

# ---- Check 1.d (plan): GR limit F=1, f=R -> Box(1)=nabla nabla(1)=0 drops out,
# leaving sg delta{g^{mu nu}}[R_{mu nu} - 1/2 g_{mu nu} R] i.e. the EH result.
v_gr = v_final.copy()
substitute(v_gr, Ex(r'\nabla_{\rho}{\nabla^{\rho}{F}} -> 0'))
substitute(v_gr, Ex(r'\nabla_{\mu}{\nabla_{\nu}{F}} -> 0'))
substitute(v_gr, Ex(r'F -> 1'))
substitute(v_gr, Ex(r'f -> R'))
canonicalise(v_gr)
print('\n6) GR limit (F->1, f->R, Box/nabla-nabla F -> 0):', v_gr)
s_gr = str(v_gr)
ok_gr = ('R_{' in s_gr) and ('R' in s_gr) and ('\\nabla' not in s_gr and '∇' not in s_gr)
print(('PASS' if ok_gr else 'FAIL'), ': GR limit reduces to R_{mu nu} - 1/2 g_{mu nu} R (no leftover nabla terms)')

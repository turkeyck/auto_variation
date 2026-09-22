"""
cadabra2 bug: free/dummy-index checker cannot see through an Accent
(\\delta{...}) boundary when summing two \\delta{...}-carrying terms in a
single Ex()-parsed string.

Bug one-line description:
    Summing two scalar terms that each carry a \\delta{...} factor inside
    one Ex()-parsed string raises "Free indices in different terms in a
    sum do not match", even when both terms are genuinely valid scalars
    (all indices properly contracted).

Expected behavior:
    `g^{\\mu\\nu} \\delta{R_{\\mu\\nu}} + \\delta{sg}` should parse without
    error -- both summands are scalars (the first has mu,nu fully
    contracted via g^{mu nu}; the second has no indices at all), so this
    is a valid sum of two scalars regardless of what sits inside the
    \\delta{...} accents.

Actual behavior:
    RuntimeError: "Free indices in different terms in a sum do not
    match." The checker appears to require the free-index pattern to
    match while still "looking through" \\delta{...} in a way that fails
    once the two terms carry structurally different content inside the
    accent (indexed vs. unindexed).

cadabra2 version tested: 2.5.14 (also confirmed 2.4.5.4 in earlier project
history, see cadabra/cadabra_utils.py). Re-verified 2026-09-12 and again
when writing this minimal reproducer (2026-09-17).

Workaround (see cadabra/cadabra_utils.py rule D1, sub_copy/replace_term):
build each additive piece as its own single-term Ex(), then combine the
already-built Ex trees with Python's `+` operator -- this bypasses the
string-level checker entirely and is demonstrated in part 2 below.

Run: wsl.exe -- python3 cadabra/bug_reports/accent_boundary_repro.py
"""
from cadabra2 import *

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu}::Indices(position=fixed).')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'R_{\mu\nu}::Symmetric.')
Ex(r'\delta{#}::Accent.')

print('--- Part 1: triggering the bug ---')
try:
    bad = Ex(r'g^{\mu\nu} \delta{R_{\mu\nu}} + \delta{sg}')
    print('NO ERROR (unexpected -- bug may be fixed in this cadabra2 version):', bad)
except RuntimeError as e:
    print('REPRODUCED:', e)

print('\n--- Part 2: the D1 workaround ---')
t1 = Ex(r'g^{\mu\nu} \delta{R_{\mu\nu}}')
t2 = Ex(r'\delta{sg}')
combined = t1 + t2
print('WORKAROUND OK, combined result:', combined)

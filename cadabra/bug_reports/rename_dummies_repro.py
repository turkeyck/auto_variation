"""
cadabra2 bug: rename_dummies() raises "No index set for index ... known"
on a valid contraction of an undeclared-property one-index tensor against
an InverseMetric/Symmetric two-index tensor, while canonicalise() handles
the identical expression without complaint.

Bug one-line description:
    `rename_dummies()` fails on `g^{mu nu} A_mu A_nu` when `A_mu` has no
    declared tensor property, even though the contraction is perfectly
    valid; `canonicalise()` on the exact same expression succeeds.

Expected behavior:
    Both `rename_dummies()` and `canonicalise()` should either succeed or
    fail consistently on the same well-formed expression -- they are
    documented as performing (at least overlapping) dummy-index-renaming
    duties.

Actual behavior:
    `rename_dummies()` raises RuntimeError: "No index set for index \\mu
    known." `canonicalise()` on the same expression succeeds silently.

cadabra2 version tested: 2.5.14 (also confirmed 2.4.5.4 in earlier project
history, see cadabra/cadabra_utils.py). Re-verified 2026-09-12 and again
when writing this minimal reproducer (2026-09-17).

Workaround (see cadabra/cadabra_utils.py rule D2): never call
rename_dummies() directly -- always use canonicalise() instead. Enforced
project-wide by cadabra/check_no_rename_dummies.py.

Run: wsl.exe -- python3 cadabra/bug_reports/rename_dummies_repro.py
"""
from cadabra2 import *

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu}::Indices(position=fixed).')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
# NOTE: A_{\mu} is deliberately left with NO declared property -- this is
# the exact condition that triggers the bug; the plain one-index tensor
# has no Symmetric/Indices-position declaration of its own.

print('--- Part 1: triggering the bug via rename_dummies() ---')
expr = Ex(r'g^{\mu\nu} A_{\mu} A_{\nu}')
print('expr before:', expr)
try:
    rename_dummies(expr)
    print('NO ERROR (unexpected -- bug may be fixed in this cadabra2 version), result:', expr)
except RuntimeError as e:
    print('REPRODUCED:', e)

print('\n--- Part 2: the D2 workaround (canonicalise() on the same expression) ---')
expr2 = Ex(r'g^{\mu\nu} A_{\mu} A_{\nu}')
try:
    canonicalise(expr2)
    print('WORKAROUND OK, result:', expr2)
except RuntimeError as e:
    print('UNEXPECTED: canonicalise() also failed:', e)

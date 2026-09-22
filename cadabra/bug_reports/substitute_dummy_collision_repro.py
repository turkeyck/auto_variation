"""
cadabra2 bug: substitute() raises a dummy-renaming error when its
replacement RHS introduces a dummy index-pair name that is already in use
as an unrelated dummy pair elsewhere in the target expression.

Bug one-line description:
    Substituting `B -> g^{mu nu} C_mu D_nu` into an expression that
    already contains `\\mu` as an unrelated dummy pair (e.g. inside
    `nabla_mu A^mu`) raises "Failed to find dummy property for $\\mu$
    while renaming dummies", even though the two `\\mu` pairs are
    logically independent (one is a new contraction introduced by the
    substitution, the other pre-existed in a completely different term).

Expected behavior:
    The substitution should succeed and produce an expression with two
    independent dummy-index pairs, however cadabra2 chooses to name them
    -- reusing a Greek letter should not be observable at the level of a
    RuntimeError.

Actual behavior:
    RuntimeError: "Failed to find dummy property for $\\mu$ while
    renaming dummies." Using a disjoint set of index letters for the
    substitution RHS (e.g. alpha, beta instead of mu, nu) avoids the
    error entirely on the identical logical substitution.

cadabra2 version tested: 2.5.14 (also confirmed 2.4.5.4 in earlier project
history, see cadabra/cadabra_utils.py). Re-verified 2026-09-12 and again
when writing this minimal reproducer (2026-09-17).

Workaround (see cadabra/cadabra_utils.py rule D3): never hand-pick index
letters for a dummy pair introduced by a substitute() RHS. Always draw
fresh names from cadabra_utils.fresh_indices(), which guarantees
distinctness from every name already handed out in the kernel session.

Run: wsl.exe -- python3 cadabra/bug_reports/substitute_dummy_collision_repro.py
"""
from cadabra2 import *

__cdbkernel__ = create_scope()

Ex(r'{\mu,\nu,\alpha,\beta}::Indices(position=fixed).')
Ex(r'g_{\mu\nu}::Metric.')
Ex(r'g^{\mu\nu}::InverseMetric.')
Ex(r'\nabla{#}::Derivative.')

print('--- Part 1: triggering the bug (reused index letters mu, nu) ---')
expr = Ex(r'B \nabla_{\mu}{A^{\mu}}')  # \mu is already a dummy pair here
print('expr before:', expr)
try:
    substitute(expr, Ex(r'B -> g^{\mu\nu} C_{\mu} D_{\nu}'))
    print('NO ERROR (unexpected -- bug may be fixed in this cadabra2 version), result:', expr)
except RuntimeError as e:
    print('REPRODUCED:', e)

print('\n--- Part 2: the D3 workaround (fresh, disjoint index letters alpha, beta) ---')
expr2 = Ex(r'B \nabla_{\mu}{A^{\mu}}')
try:
    substitute(expr2, Ex(r'B -> g^{\alpha\beta} C_{\alpha} D_{\beta}'))
    print('WORKAROUND OK, result:', expr2)
except RuntimeError as e:
    print('UNEXPECTED: fresh-index substitution also failed:', e)
